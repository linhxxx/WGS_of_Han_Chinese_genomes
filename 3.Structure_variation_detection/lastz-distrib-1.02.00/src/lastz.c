//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: lastz.c
//
//----------

char* programName            = "lastz";
char* programVersionMajor    = VERSION_MAJOR;
char* programVersionMinor    = VERSION_MINOR;
char* programVersionSubMinor = VERSION_SUBMINOR;
char* programRevisionDate    = REVISION_DATE;

char* svnRevisionNumber      = SUBVERSION_REV;

//----------
//
// lastz-- Local Alignment Search Tool, blastZ-like
//	Find pairwise local alignments between a target DNA sequence and a series
//	of query sequences.  Query sequences can be DNA or quantum DNA.
//
// GOAL 1: Align sequences as large as 250M bases (about the size of human
//         chromsome 1) on a workstation with 32-bit addressing and 2G byte of
//         memory, with a seed up to weight 24 (equivalent to a 12-mer).
//
// GOAL 2: It should be relatively easy to try out different alignment
//         strategies.
//
//----------
//
// Algorithmic Overview
//
// Caveat:  many details are left out or glossed over in this description.  In
// some cases accuracy has been sacrificed for clarity.
//
// The algorithm consists of the following stages:
//	1 (SEED)   discovery of short near-matches ("seed hits")
//	2 (HSP)    extension to ungapped High-scoring Segment Pairs (HSPs)
//	3 (CHAIN)  reduction of HSPs to the best-scoring "chain"
//	4 (ANCHOR) reduction of (remaining) HSPs to single positions ("anchors")
//	5 (ALIGN)  extension of anchors to gapped alignments
//	6 (INTERP) "interpolation" between alignments at a higher sensitivity
//
// Most stages are optional.  For example, you can leave out chaining, or you
// can stop as soon as HSPs are found.  However, there are some dependencies
// between stages.  SEED and HSP are performed together and are required for
// anything else, ANCHOR and ALIGN are performed together, and INTERP requires
// ALIGN.  INTERP can be thought of as repeating the first five stages on
// whatever failed to align.
//
// The process is repeated, (mostly) independently, for every query sequence.
// The exceptions are a seed position table is built once, and dynamic masking
// has interactions between queries. (These are described below).
//
// The SEED stage makes use of a large table containing the positions (in
// the target sequence) of every possible seed bit pattern.  For this discussion
// seed can be thought of as 12-mers, but we support more general seeds (see
// seeds.c).  The table is indexed by the seed bit pattern, which for an 12-mer
// is a 24 bit value.  For each seed value the table gives a list of all
// positions where that seed value can be found.  The table requires 4*(L+4^W)
// bytes, where W is the seed length and L is the sequence length.  For W=12
// and L = 250MB, this is a little larger than 1GB.  For details on how the
// table is stored, see pos_table.c.
//
// Having created the position table, the SEED stage scans a query sequence and
// looks up matches ("seed hits") for every query seed.  As each seed hit is
// found, the HSP stage is performed, extending it along the diagonal in both
// directions until the score drops off.  HSPs that do not meet the score
// threshold are discarded (strictly speaking, they are not HSPs).  To improve
// performance, we keep track of how far we have progressed along each diagonal,
// and quickly discard any seed hits that fall into previously identified HSPs.
// Since a complete diagonal tracking array could be huge (for two 250M base
// sequences it would be 4 * 500M = 2G), we use a much smaller array and hash
// diagonals into its index space.  This results in a small loss of sensitivity,
// as hash collisions between diagonals can cause seed hits on other diagonals
// to be hidden by HSPs on other diagonals.  See seed_search.c and diag_hash.c.
//
// The CHAIN stage finds the highest scoring series of HSPs in which each HSP
// begins strictly before the start of the next.  All HSPs not on this chain
// are discarded.  This is useful when the query and target are known to be
// syntenic.  The algorithm processes the HSPs in order along the target, build
// chains by adding the next HSP to the best previous viable chain.  See
// chain.c.
//
// The ANCHOR stage reduces each HSP to a single point.  A constant-width
// window is slid across the HSP and the midpoint of the highest-scoring window
// is chosen as the anchor.
//
// The ALIGN stage extends anchors into gapped alignments.  One-sided extension
// is performed in two directions from the anchor point, the two resulting
// alignments are joined at the anchor, and if the score is high enough, this
// becomes an alignment in the output file.  Extension is computed per the
// standard (but optimized) 3-state affine gap dynamic programming recurrence.
// The DP matrix is evaluated only in a small region straddling the high-scoring
// path, with only enough memory to store the widest row of that region.  See
// gapped_extend.c.
//
// Dynamic masking is (optionally) performed during the ALIGN stage.  We keep a
// count of how many times each target base has been in an alignment.  When a
// base reaches a threshold, it is assumed to be a repeat and is masked from the
// sequence (by chaing it to an 'x').  This only affects subsequent queries in
// the same run.  See masking.c.
//
// The INTERP stage repeats all of the previous steps in the regions between
// those gapped alignments.  As of this writing the interpolation seed is a
// 7-mer exact match, for compatibility with BLASTZ.
//
// Input formats can be FASTA, nib, or qdna files, but other formats can be
// added without too much pain.  See sequences.c.
//
// A variety of pairwise output formats are supported.  The primary format is
// LAV for compatibility with BLASTZ.  See lav.c, gfa.c, axt.c, maf.c, cigar.c,
// genpaf.c, text_align.c and align_diffs.c.  It is relatively easy to add other
// formats (although it would be even easier if this were written in an object-
// oriented language).
//
//----------

//----------
//
// other files
//
//----------

#include <stdlib.h>				// standard C stuff
#define  true  1
#define  false 0
#include <stdio.h>				// standard C i/o stuff
#include <string.h>				// standard C string stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include <math.h>				// standard C math stuff
#include <time.h>				// standard C time stuff
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "seeds.h"				// seed strategy stuff
#include "pos_table.h"			// position table stuff
#include "capsule.h"			// multi-process sharing stuff
#include "seed_search.h"		// seed hit search stuff
#include "quantum.h"			// quantum DNA search stuff
#include "segment.h"			// segment table management stuff
#include "chain.h"				// segment chaining stuff
#include "gapped_extend.h"		// gapped alignment stuff
#include "tweener.h"			// interpolated alignment stuff
#include "masking.h"			// dynamic masking stuff
#include "infer_scores.h"		// scoring inference stuff
#include "edit_script.h"		// alignment edit script stuff
#include "diag_hash.h"			// diagonals hashing stuff
#include "identity_dist.h"		// identity distribution stuff
#include "coverage_dist.h"		// query coverage distribution stuff
#include "continuity_dist.h"	// query continuity distribution stuff
#include "output.h"				// alignment outout format stuff
#include "maf.h"				// maf alignment format stuff
#include "genpaf.h"				// genpaf alignment format stuff
#include "text_align.h"			// textual alignment format stuff

#define  lastz_owner			// (make this owner of program-wide globals)
#include "lastz.h"				// lastz program-wide stuff

#define helpout stdout			// stream to write help messages to

//----------
//
// stuff for crude profiling
//
//----------

//#define useStandardClock  (define at build time)

#ifdef useStandardClock
#define read_clock() clock()
#define clocksPerSec CLOCKS_PER_SEC
#endif // useStandardClock

#ifndef useStandardClock
#include <sys/time.h>
#define read_clock() microsec_clock()
#define clocksPerSec 1000000
u64 microsec_clock (void);
u64 microsec_clock (void)
	{
	static int		failed = false;
	struct timeval	time1;
	int				err;

	if (failed)	return 0;	// (previous call to gettimeofday has failed)
	err = gettimeofday (&time1, NULL);
	if (err != 0) { failed = true;  return 0; }

	return (((u64) time1.tv_sec) * 1000000) + time1.tv_usec;
	}
#endif // not useStandardClock

// clock for dbgQueryProgress

s64 dbgQueryProgressClock = 0;

// build-specific profiling stuff

#ifndef dbgTiming
#define dbg_timing_sub(v) ;
#define dbg_timing_add(v) ;
#define dbg_timing_copy(dst,src) ;
#define dbg_timing_report(v,s) ;
#endif // not dbgTiming

#ifdef dbgTiming

s64 debugClockTotal = 0,
    debugClockSeq1 = 0,
    debugClockSeq2 = 0,
    debugClockPosTable = 0,
    debugClockQueryTotal = 0,
    debugClockSegTable = 0,
    debugClockChaining = 0,
    debugClockGappedExtend = 0,
    debugClockInterpolation = 0,
    debugClockOutput = 0;

#define dbg_timing_sub(v)  { v -= (s64) read_clock();  }
#define dbg_timing_add(v)  { v += (s64) read_clock();  }
#define dbg_timing_copy(dst,src) { dst = src; }

#define dbg_timing_report(v,s) { fprintf(stderr,"%-26s %.3f\n",s":",((float)(v))/clocksPerSec); }
#endif // dbgTiming

//----------
//
// private global data
//
//----------

static int showProgress = false;

// we keep two copies of the control data, one for the primary alignment, and
// another for alignments used to infer a scoring set;  for historical reasons
// we call these "lz" (short for lastz) and "iz" (short for inferz);  we keep a
// pointer to whichever of these is 'active' at any given time

static control  lzParams;
static control  izParams;

control* currParams;	// (nota bene:  currParams is accessed by output.c)

// command line options

// defaults

static const control defaultParams =
	{
	NULL,NULL,							// seq1, seq1Filename
	NULL,NULL,							// seq2, seq2Filename
	NULL,NULL,							// rev1, rev2

	false,false,						// inferScores, inferOnly
	{
	NULL,								// ic.inferFilename
	100,true,							// ic.inferScale, ic.writeAsInt
	ratioNone,							// ic.hspThresholdIsRatio,
	ratioNone,							// ic.gappedThresholdIsRatio
	ratioNone,							// ic.gapOpenIsRatio
	ratioNone,							// ic.gapExtendIsRatio
	30,0,								// ic.subIterations, ic.gapIterations
	false,								// ic.idIsPercentile
	},
	false,								// selfCompare

	true,								// doSeedSearch
	(s8*) nuc_to_bits,					// charToBits
	(s8*) upper_nuc_to_bits,			// upperCharToBits
	1,									// whichStrand
	1,									// step

	NULL,								// hitSeed (default is defaultSeedString)
	28,									// maxIndexBits
	1,									// withTrans
	false,								// noHitFiltering
	0,0,								// twinMinSpan, twinMaxSpan (trumped
										// .. by defaultTwinsYes, etc. below)
	hitSimple,							// basicHitType
	-1,-1,								// minMatches, maxTransversions
	false,								// filterCaresOnly
#ifndef noSeedHitQueue
	defaultSeedHitQueueSize,			// seedHitQueueSize
#endif // not noSeedHitQueue

	false,false,						// readCapsule, writeCapsule,
	NULL,NULL,NULL,						// capsuleFile, capsuleFilename, capsule

	NULL,NULL,							// anchorsFile, anchorsFilename

	gfexXDrop,							// gfExtend
	false,								// mergeAnchors
	false,0,0,							// chain, chainDiag, chainAnti
	true,								// gappedExtend

	NULL,NULL,							// scoring, maskedScoring
	0,									// xDrop
	0,									// yDrop
	false,false,						// xDropUntrimmed, yDropUntrimmed
	{'S',3000,0,0},						// hspThreshold
	{'S',0,   0,0},						// gappedThreshold
	true,false,							// entropicHsp, reportEntropy
	false,								// gappedAllBounds
	-1,									// mirror (not used, default hardcoded)
	false,								// inhibitTrivial
	80*1024*1024,						// tracebackMem
	NULL,								// traceback
	false,false,						// nIsAmbiguous, allowAmbiDNA
	false,0,							// hspToGapped,searchLimit

	0.0,0,								// wordCountKeep, wordCountLimit
	0,									// maxWordCountChasm
	0,									// dynamicMasking
	false,								// reportCensus
	NULL,								// censusFile
	NULL,								// censusFilename
	'B',								// censusKind

	0.0,1.0,							// minIdentity, maxIdentity
	0.0,1.0,							// minCoverage, maxCoverage
	0.0,1.0,							// minContinuity, maxContinuity
	0,									// minMatchCount
#ifdef densityFiltering
	0.0,								// maxDensity
#endif // densityFiltering
	NULL,NULL,							// outputFilename, outputFile
	fmtLav,NULL,						// outputFormat, outputInfo
	false,								// endComment
	false,								// needTrueLengths
	false,								// deGapifyOutput
	NULL,NULL,							// dotplotFilename, dotplotFile

	0,									// innerThreshold
	NULL,								// innerSeed
	20000,								// innerWindow

	false,false,						// targetIsQuantum, queryIsQuantum
	-1,									// ballScore

	true,								// lajCompatible
	0,									// textContext
	NULL,								// args
	0,									// verbosity
	false,								// reportTiming
	false,								// reportStats
	false,								// showStats
	NULL,								// statsFile
	NULL,								// statsFilename
	spt_dont							// showPosTable
	};

static const char* defaultSeedString = seed_12of19;

static const int   defaultTwinsYes   = false;
static const int   defaultTwinMinGap = 0;
static const int   defaultTwinMaxGap = 10;

static int dbgShowMatrix    = false;
static int dbgShowParams    = false;
static int dbgShowHsps      = false;
static int dbgShowAnchors   = false;
static int dbgSortAnchorsByDiag = false;
static int dbgInhibitSegmentReduction = false;
static int dbgMasking       = false;
static int dbgQueryProgress = 0;

#define innerWordSize	7				// word size for inner alignment seed
										// .. hits

static const float defaultBallScoreFactor = 0.75;

// anchors--
//	Whenever the process will go beyond just finding gap-free extensions, the
//	segments that will become anchors (e.g. HSPs) are collected in this table
//	instead of being written to the console.

#define numDefaultAnchors 4000

static segtable* anchors      = NULL;
static segtable* savedAnchors = NULL;

// miscellany

#define chainScale 100

sthresh scratchThreshold = {'S',0,0,0};

//----------
//
// prototypes for private functions
//
//----------

int main (int argc, char** argv);

static seq*  capsule_target          (capinfo* cap, u8** targetRev);
static postable* capsule_position_table (capinfo* cap, seq* seq,
                                      seed* hitSeed, u32 step);
static score chain_connect_penalty   (segment* seg1, segment* seg2, int scale);
static void  remove_interval_seeds   (unspos b, unspos e, void* info);
static u32   report_hsps             (void* info,
                                      unspos pos1, unspos pos2, unspos length,
                                      score s);
static u32   collect_hsps            (void* info,
                                      unspos pos1, unspos pos2, unspos length,
                                      score s);
static void  parse_options           (int argc, char** argv,
                                      control* lzParams, control* izParams);
static void  create_seed_structure   (control* lzParams, char** seedString,
                                      int haveWithTrans, int twinsYes,
                                      int minGap, int maxGap);

static void  print_params            (FILE* f, control* lzParams);

static void  read_control_file_by_name (char* name, control* params);
static void  read_control_file       (FILE* f, char* name, control* params);

static void  print_options           (void);

static int   name_spec_is_quantum    (char* spec);

static void  lastz_zero_stats        (void);
static void  lastz_show_stats_before (FILE* f);
static void  lastz_show_stats        (FILE* f);

//----------
//
// lastz--
//  Main program
//
//----------

int main
   (int				argc,
	char**			argv)
	{
	FILE*			statsF        = NULL;
	seq*			target        = NULL;
	seq*			query         = NULL;
	postable*		targPositions = NULL;
	int				freeTargetRev = false;
	u8*				targetRev     = NULL;
	tback*			traceback     = NULL;
	census*			targCensus    = NULL;
	hitprocessor	hitProc;
	void*			hitProcInfo;
	time_t			startClock, endClock;
	unspos			coverageLimit;
	int				reverseNeeded, hspsAreAdaptive, tableWillBeUsed, queryExists;
	int				numQueries = 0;
	u8				rCh, cCh;

	dbg_timing_sub (debugClockTotal);
	dbgQueryProgressClock = -((s64) read_clock());

	debug = 0;

	lastz_zero_stats         ();
	pos_table_zero_stats     ();
	capsule_zero_stats       ();
	seed_search_zero_stats   ();
	quantum_zero_stats       ();
	chain_zero_stats         ();
	gapped_extend_zero_stats ();
	tweener_zero_stats       ();
	masking_zero_stats       ();
	infer_scores_zero_stats  ();

	//////////
	// fetch arguments
	//////////

	currParams = NULL;
	parse_options (argc, argv, &lzParams, &izParams);
	currParams = &lzParams;

	// open stats file

	if (currParams->statsFilename != NULL)
		{
		currParams->statsFile = fopen_or_die (currParams->statsFilename, "wt");
		free_if_valid ("stats file name", currParams->statsFilename);
		}

	if (currParams->showStats)
		statsF = (currParams->statsFile != NULL)? currParams->statsFile : stderr;

	if ((currParams->inferScores) && (statsF != NULL))
		izParams.statsFile = statsF;

	// open anchors file;  note that we don't dispose of the file name here, so
	// that we can use it in error messages later

	if (currParams->anchorsFilename != NULL)
		{
		currParams->anchorsFile  = fopen_or_die (currParams->anchorsFilename, "rt");
		currParams->mergeAnchors = true;
		}

	//////////
	// open the sequence files
	//
	// We also load the first sequence from the input files here.  This allows
	// us to check and report problems (such as an extra target or an empty
	// query) immediately, rather than wasting what might turn out to be a lot
	// of processing time.
	//////////

	startClock = clock();

	// open and load target

	memory_checkpoint ("[[* Opening Files ]]\n");

	dbg_timing_sub (debugClockSeq1);

	reverseNeeded = ((currParams->gappedExtend)
	              || ((currParams->inferScores) && (izParams.gappedExtend)));

	if (lzParams.capsule != NULL)
		{
		if (!reverseNeeded)
			target = capsule_target (lzParams.capsule, NULL);
		else
			{
			target = capsule_target (lzParams.capsule, &targetRev);
			freeTargetRev = false;
			currParams->rev1 = targetRev;
			}
		currParams->seq1 = target;
		}
	else
		{
		if (currParams->verbosity >= 5)
			fprintf (stderr, "opening \"%s\"\n", currParams->seq1Filename);

		target = open_sequence_file (currParams->seq1Filename, seq_type_unknown,
		                             currParams->needTrueLengths,
		                             currParams->allowAmbiDNA, NULL);

		if ((currParams->targetIsQuantum) && (target->fileType != seq_type_qdna))
			suicidef ("%s does not contain quantum DNA", target->fileName);

		currParams->seq1 = target;

		if ((!load_sequence (target)) || (target->len == 0))
			suicidef ("%s contains no sequence", target->fileName);
		if (another_sequence (target))
			suicidef ("%s contains more than one sequence\n"
			          "consider using the \"multiple\" action (see \"lastz --help=files\")",
			          target->fileName);

		if (reverseNeeded)
			{
			targetRev = (u8*) copy_reverse_of_string ((char*) target->v, target->len);
			freeTargetRev = true;
			currParams->rev1 = targetRev;
			}
		}

	dbg_timing_add (debugClockSeq1);

	if ((currParams->anchorsFile != NULL) && (target->revCompFlags != rcf_forward))
		suicidef ("can't use --segments with reverse-complement of target (%s)",
		          target->fileName);

	if ((currParams->dynamicMasking > 0) || (currParams->reportCensus))
		targCensus = new_census (target->len, currParams->censusKind, currParams->dynamicMasking);

	// now that we know the target length, set the hsp and gapped thresholds
	// if they are a percentage of the target length

	resolve_score_thresh (&currParams->hspThreshold,    target->len);
	resolve_score_thresh (&currParams->gappedThreshold, target->len);

	if (currParams->inferScores)
		{
		resolve_score_thresh (&izParams.hspThreshold,    target->len);
		resolve_score_thresh (&izParams.gappedThreshold, target->len);
		}

	// open and load query

	if ((currParams->doSeedSearch) || (currParams->inferScores))
		{
		if (currParams->verbosity >= 5)
			{
			if (currParams->seq2Filename != NULL)
				fprintf (stderr, "opening \"%s\"\n", currParams->seq2Filename);
			else
				fprintf (stderr, "opening unnamed query file\n");
			}

		dbg_timing_sub (debugClockSeq2);
		if (currParams->inferScores)
			query = open_rewindable_sequence_file (currParams->seq2Filename,
			                                       seq_type_unknown,
			                                       currParams->needTrueLengths,
			                                       currParams->allowAmbiDNA,
			                                       NULL);
		else
			query = open_sequence_file            (currParams->seq2Filename,
			                                       seq_type_unknown,
			                                       currParams->needTrueLengths,
			                                       currParams->allowAmbiDNA,
			                                       currParams->maskedScoring->qToComplement);

		if ((currParams->queryIsQuantum) && (query->fileType != seq_type_qdna))
			suicidef ("%s does not contain quantum DNA", query->fileName);

		currParams->seq2 = query;
		if (!another_sequence (query))
			suicidef ("%s contains no sequence", query->fileName);

		dbg_timing_add (debugClockSeq2);
		}

	if (dbgShowParams)
		print_params (stderr, &lzParams);

	// check for bad combination of partitioned sequence vs other options

	if ((target->partition.p != NULL)
	 || ((query != NULL) && (query->partition.p != NULL)))
		{
		if ((currParams->doSeedSearch)
		 && (!currParams->inferOnly)
		 && ((currParams->outputFormat == fmtGfa)
		  || (currParams->outputFormat == fmtGfaNoScore)))
			suicide ("multiple action cannot be used with --gfa");

		if ((currParams->doSeedSearch)
		 && (!currParams->inferOnly)
		 && ((currParams->outputFormat == fmtLav)
		  || (currParams->outputFormat == fmtLavComment)
		  || (currParams->outputFormat == fmtLavScore)
		  || (currParams->outputFormat == fmtLavText)))
			suicide ("multiple action cannot be used with --lav\n"
			         "(lav has requirements on the order of alignments that would require additional\n"
			         " computation;  use \"--help=formats\" to see other options for output)");
		}

	if ((target->partition.p != NULL) && (currParams->selfCompare))
		suicide ("multi cannot be used with --self");

	// allocate traceback memory

	if ((currParams->gappedExtend)
	 || ((currParams->inferScores) && (izParams.gappedExtend)))
		{
		traceback = new_traceback (currParams->tracebackMem);
		currParams->traceback = traceback;
		}

	//////////
	// build a position table for the target sequence
	//
	// If we are inferring scores, we'd like to share the same position table
	// for both inference and for the final alignment.  However, it is natural
	// to use a less sensitive seed or step during inference.  So in case the
	// seed or step is different, we use the inference control values here and
	// will rebuild the table after inference is finished, if needed.
	//
	// If we are to get the position table from a capsule file, we do so here.
	// Note that scoring inference and reading a capsule file are mutually
	// exclusive options.
	//////////

	if (currParams->anchorsFile == NULL)
		{
		dbg_timing_sub (debugClockPosTable);

		if (currParams->verbosity >= 1)
			fprintf (stderr, "building position table for %s\n",
			                 target->fileName);

		if (lzParams.capsule != NULL)
			targPositions = capsule_position_table
			                   (lzParams.capsule, target,
			                    currParams->hitSeed, currParams->step);
		else if (currParams->inferScores)
			targPositions = build_seed_position_table
							   (target, 0, target->len,
								currParams->upperCharToBits, izParams.hitSeed,
								izParams.step);
		else if (!currParams->targetIsQuantum)
			{
			targPositions = build_seed_position_table
							   (target, 0, target->len,
								currParams->upperCharToBits, currParams->hitSeed,
								currParams->step);
			if (currParams->wordCountKeep > 0)
				currParams->wordCountLimit = find_position_table_limit
				                                (targPositions,
				                                 currParams->wordCountKeep);
			if (currParams->wordCountLimit > 0)
				limit_position_table (targPositions, currParams->wordCountLimit,
				                      currParams->maxWordCountChasm);
			}
		else // if (currParams->targetIsQuantum)
			{
			targPositions = build_quantum_seed_position_table
							   (target, 0, target->len,
							    currParams->maskedScoring->bottleneck,
							    currParams->maskedScoring->qToBest,
								currParams->hitSeed, currParams->step);
			if (currParams->wordCountKeep > 0)
				currParams->wordCountLimit = find_position_table_limit
				                                (targPositions,
				                                 currParams->wordCountKeep);
			if (currParams->wordCountLimit > 0)
				limit_position_table (targPositions, currParams->wordCountLimit,
				                      currParams->maxWordCountChasm);
			}

		dbg_timing_add (debugClockPosTable);

		lastz_show_stats_before (statsF);
		pos_table_show_stats (statsF, targPositions);
		}

	//////////
	// perform scoring inference, if requested
	//////////

	if (currParams->inferScores)
		{
		// switch control to the inference parameters

		currParams = &izParams;
		izParams.seq1Filename = lzParams.seq1Filename;
		izParams.seq1         = lzParams.seq1;
		izParams.seq2Filename = lzParams.seq2Filename;
		izParams.seq2         = lzParams.seq2;

		// perform the scoring inference (note that this will rewind the query
		// file)

		if (currParams->hspThreshold.t == 'S') coverageLimit = 0;
		else coverageLimit = currParams->hspThreshold.c;

		anchors = new_segment_table (numDefaultAnchors, coverageLimit);

		lzParams.scoring = drive_scoring_inference
		                       (currParams,
		                        target, targetRev, targPositions, query,
		                        traceback);

		// switch control back to the main parameters

		currParams = &lzParams;

		// fill in score-based parameters

		currParams->maskedScoring = masked_score_set (currParams->scoring);

		if (currParams->xDrop < 0)
			{
			rCh = currParams->scoring->rowChars[0];
			cCh = currParams->scoring->colChars[0];
			currParams->xDrop = 10 * currParams->scoring->sub[rCh][cCh];
			}

		if (currParams->yDrop < 0)
			currParams->yDrop = currParams->scoring->gapOpen + 300 * currParams->scoring->gapExtend;

		// rebuild the target sequence position table if it would be different

		tableWillBeUsed = ((currParams->doSeedSearch)
		                || (currParams->showPosTable != spt_dont)
		                || (currParams->writeCapsule));

		if ((tableWillBeUsed)
		 && (!is_same_seed (izParams.hitSeed, currParams->hitSeed))
		 && (izParams.step != currParams->step))
			{
			free_position_table (targPositions);
			if (!currParams->targetIsQuantum)
				targPositions = build_seed_position_table
								   (target, 0, target->len,
									currParams->upperCharToBits, currParams->hitSeed,
									currParams->step);
			else // if (currParams->targetIsQuantum)
				targPositions = build_quantum_seed_position_table
								   (target, 0, target->len,
								    currParams->maskedScoring->bottleneck,
								    currParams->maskedScoring->qToBest,
									currParams->hitSeed, currParams->step);
			if (currParams->wordCountKeep > 0)
				currParams->wordCountLimit = find_position_table_limit
				                                (targPositions,
				                                 currParams->wordCountKeep);
			if (currParams->wordCountLimit > 0)
				limit_position_table (targPositions, currParams->wordCountLimit,
				                      currParams->maxWordCountChasm);
			}
		}

	if (currParams->showPosTable == spt_distribution)
		{
		poscount* posDist = position_table_count_distribution (targPositions);
		poscount* pd;
		fprintf (currParams->outputFile,
	            "seed-word counts distribution table for %s:\n",
		        currParams->seq1->fileName);
		for (pd=posDist ; pd->occurrences!=0 ; pd++)
			fprintf (currParams->outputFile, unsposFmt " " unsposFmt "\n",
			         pd->count, pd->occurrences);
		free_if_valid ("seed word position counts distribution",  posDist);
		}
	else if (currParams->showPosTable != spt_dont)
		{
		if (currParams->showPosTable == spt_table)
			fprintf (currParams->outputFile,
		            "seed-word positions table for %s:\n",
			        currParams->seq1->fileName);
		else if (currParams->showPosTable == spt_countsonly)
			fprintf (currParams->outputFile,
		            "seed-word counts table for %s:\n",
			        currParams->seq1->fileName);
		else // if (currParams->showPosTable == spt_withcounts)
			fprintf (currParams->outputFile,
		            "seed-word counts and positions table for %s:\n",
			        currParams->seq1->fileName);
		dump_position_table (currParams->outputFile,
		                     targPositions, currParams->hitSeed,
		                        (currParams->showPosTable == spt_table)
		                     || (currParams->showPosTable == spt_withcounts),
		                        (currParams->showPosTable == spt_countsonly)
		                     || (currParams->showPosTable == spt_withcounts));
		printf ("\n");
		}

	//////////
	// if we are only writing a capsule file, do so and quit
	//////////

	if (currParams->writeCapsule)
		{
		u64 capSize;

		currParams->capsuleFile = fopen_or_die (currParams->capsuleFilename, "wb");
		capSize = write_capsule_file (currParams->capsuleFile,
									  currParams->capsuleFilename,
		                              target, targetRev, targPositions,
		                              currParams->hitSeed);
		fclose_if_valid (currParams->capsuleFile);
		currParams->capsuleFile = NULL;
		endClock = clock();
		printf ("%s byte target sequence capsule written to %s\n",
		        unitize(capSize), currParams->capsuleFilename);
		goto show_stats_and_clean_up;
		}

	//////////
	// perform the alignment of the query (or queries) to the target
	//////////

	dbg_timing_copy (debugClockQueryTotal, debugClockSeq2);
	dbg_timing_sub  (debugClockQueryTotal);
	numQueries = 0;

	// if the user doesn't want the alignment we can quit now

	if (((!currParams->inferScores) && (!currParams->doSeedSearch))
	 || (( currParams->inferScores) && (currParams->inferOnly)))
		{
		endClock = clock();
		goto show_stats_and_clean_up;
		}

	// initialize (or empty) the anchors list

	if (currParams->hspThreshold.t == 'S') coverageLimit = 0;
	else coverageLimit = currParams->hspThreshold.c;

	if (anchors == NULL)
		anchors = new_segment_table (numDefaultAnchors, coverageLimit);
	else
		{
		empty_segment_table (anchors);
		limit_segment_table (anchors, coverageLimit);
		}

	if (currParams->anchorsFile == NULL)
		set_up_hit_processor (currParams, (targCensus!=NULL),
		                      &hitProc, &hitProcInfo);
	else
		{
		hitProc     = NULL;
		hitProcInfo = NULL;
		}

	// search for hits in each query sequence

	hspsAreAdaptive = (currParams->hspThreshold.t != 'S');

	print_job_header ();
	print_options ();

	while (true)
		{
		memory_checkpoint_1 ("[[* Query #%d (loading)]]\n", numQueries);

		dbg_timing_sub (debugClockSeq2);
		queryExists = load_sequence (query);
		dbg_timing_add (debugClockSeq2);
		if (!queryExists) break;
		numQueries++;

		if (query->len == 0) continue;

		memory_checkpoint_2 ("[[* Query #%d, %s (loaded) ]]\n",
		                     numQueries,
		                     (query->useFullNames)? query->header
		                                          : query->shortHeader);

		if (currParams->whichStrand < 0)
			rev_comp_sequence (query, currParams->scoring->qToComplement);

		if ((dbgQueryProgress != 0)
		 && ((dbgQueryProgress == 1) || (numQueries % dbgQueryProgress == 1)))
		 	{
			dbgQueryProgressClock += (s64) read_clock();
			fprintf (stderr, "(%.3fs) ", ((float)(dbgQueryProgressClock))/clocksPerSec);
			dbgQueryProgressClock = -((s64) read_clock());

			fprintf (stderr, "processing query %d", numQueries);
			if (query->shortHeader != NULL)
				fprintf (stderr, ": %s", query->shortHeader);
			else if (query->header != NULL)
				fprintf (stderr, ": %s", query->header);
			fprintf (stderr, "\n");
			}

		init_output_for_query ();

		// search for "forward" strand hits (and finish alignment for them if
		// hsps are non-adaptive)

		if (currParams->verbosity >= 1)
			fprintf (stderr, "searching for matches in %s%s\n",
							 query->fileName,
							 (currParams->whichStrand>=0)? ""
		                                             : ", (reverse strand)");

		if (currParams->hspToGapped)
			{
			hitprocsimple* simpleInfo  = (hitprocsimple*) hitProcInfo;
			hitrepgappily* gappilyInfo = (hitrepgappily*) simpleInfo->hp.reporterInfo;

			if (currParams->rev2 != NULL)
				suicidef ("internal error, currParams->rev2 is not NULL");

			currParams->rev2 = (u8*) copy_reverse_of_string ((char*) query->v, query->len);

			gappilyInfo->seq2 = currParams->seq2;
			gappilyInfo->rev2 = currParams->rev2;
			}

		start_one_strand (target, targPositions, query, true,
		                  hitProc, hitProcInfo);
		if (!hspsAreAdaptive)
			finish_one_strand (target, targetRev, targPositions, query, NULL,
							   traceback, targCensus);

		// search for reverse strand hits

		if (currParams->whichStrand > 0)
			{
			if (currParams->verbosity >= 1)
				fprintf (stderr, "searching for matches in %s, (reverse strand)\n",
								 query->fileName);

			rev_comp_sequence (query, currParams->scoring->qToComplement);

			if (currParams->hspToGapped)
				strncpy_reverse (/* to   */ (char*) currParams->rev2,
				                 /* from */ (char*) query->v,
				                 /* size */         query->len);

			start_one_strand  (target, targPositions, query, !hspsAreAdaptive,
			                   hitProc, hitProcInfo);

			// finish alignment for reverse strand hits

			if (hspsAreAdaptive) split_anchors (query->revCompFlags);
			finish_one_strand (target, targetRev, targPositions, query, NULL,
							   traceback, targCensus);


			if (hspsAreAdaptive) swap_anchors ();
			if (hspsAreAdaptive)
				rev_comp_sequence (query, currParams->scoring->qToComplement);
			}

		// finish alignment for "forward" strand hits

		if (hspsAreAdaptive)
			finish_one_strand (target, targetRev, targPositions, query, NULL,
							   traceback, targCensus);

		if ((currParams->hspToGapped) && (currParams->rev2 != NULL))
			{
			free_if_valid ("reverse query (currParams->rev2)", currParams->rev2);
			currParams->rev2 = NULL;
			}
		}

	if (currParams->anchorsFile != NULL)
		{
		// make sure all segments were read

		read_segment_table (currParams->anchorsFile, currParams->anchorsFilename,
		                    NULL, NULL, NULL);
		}

	endClock = clock();

	print_m_stanza (targCensus);
	if (currParams->reportCensus)
		{
		u32 savedThresh = targCensus->maskThresh;
		targCensus->maskThresh = 0;
		if (currParams->censusFilename == NULL)
			print_census_stanza (targCensus);
		else
			{
			if (currParams->censusFile == NULL)
				currParams->censusFile = fopen_or_die (currParams->censusFilename, "wt");
			print_census (currParams->censusFile, target, targCensus, '\t');
			}
		targCensus->maskThresh = savedThresh;
		}

#ifdef collect_stats
	if (currParams->gfExtend != gfexNoExtend)
		{
		print_generic (currParams->outputFile,
		               "gap_free_extensions=%jd",
		               seed_search_hsps()
		             + seed_search_low_scoring_hsps());
		if (!currParams->reportStats)
			print_generic (currParams->outputFile,
			               "bp_extended=%jd",
			               seed_search_bp_extended());
		}
#endif // collect_stats

	if (currParams->reportStats)
		{
		seed_search_generic_stats   (currParams->outputFile, print_generic);
		if (query->fileType == seq_type_qdna)
			quantum_generic_stats   (currParams->outputFile, print_generic);
		chain_generic_stats         (currParams->outputFile, print_generic);
		gapped_extend_generic_stats (currParams->outputFile, print_generic);
		if (currParams->innerSeed != NULL)
			tweener_generic_stats   (currParams->outputFile, print_generic);
		if (currParams->dynamicMasking > 0)
			masking_generic_stats (currParams->outputFile, print_generic);
		if (currParams->inferScores)
			infer_scores_generic_stats (currParams->outputFile, print_generic);
		}

	print_job_footer ();

show_stats_and_clean_up:
	{
	float runTime = ((float)(endClock-startClock))/CLOCKS_PER_SEC;
	if (currParams->reportTiming)
		print_generic (currParams->outputFile, "runtime=%.3f", runTime);
	lastz_set_stat (runTime, runTime);
	}

	capsule_show_stats       (statsF);
	seed_search_show_stats   (statsF);
	if ((query != NULL) && (query->fileType == seq_type_qdna))
		quantum_show_stats   (statsF);
	chain_show_stats         (statsF);
	gapped_extend_show_stats (statsF);
	if (currParams->innerSeed != NULL)
		tweener_show_stats (statsF);
	if (currParams->dynamicMasking > 0)
		{
		pos_table_show_stats_after (statsF);
		masking_show_stats (statsF);
		}
	if (currParams->inferScores)
		infer_scores_show_stats (statsF);
	lastz_show_stats (statsF);

	if (currParams->endComment)
		print_end_comment ();

	//////////
	// clean up
	//
	// note that we don't bother to dispose of allocated memory unless we are
	// going to be running the valgrind memory checker
	//////////

	memory_checkpoint ("[[* Cleanup ]]\n");


#if ((defined trackMemoryUsage) || (defined valgrindMemoryCheck))

	free_if_valid        ("lz.outputFilename",  lzParams.outputFilename);   lzParams.outputFilename   = NULL;
	fclose_if_valid      (lzParams.outputFile);                             lzParams.outputFile       = NULL;
	free_if_valid        ("lz.dotplotFilename", lzParams.dotplotFilename);  lzParams.dotplotFilename  = NULL;
	fclose_if_valid      (lzParams.dotplotFile);                            lzParams.dotplotFile      = NULL;
	free_if_valid        ("lz.seq1Filename",    lzParams.seq1Filename);     lzParams.seq1Filename     = NULL;
	free_sequence        (lzParams.seq1);                                   lzParams.seq1             = NULL;
	if (freeTargetRev)
		{ free_if_valid  ("targetRev",          targetRev);                 targetRev                 = NULL; }
	free_if_valid        ("lz.seq2Filename",    lzParams.seq2Filename);     lzParams.seq2Filename     = NULL;
	free_sequence        (lzParams.seq2);                                   lzParams.seq2             = NULL;
	free_if_valid        ("lz.args",            lzParams.args);             lzParams.args             = NULL;
	free_score_set       ("lz.scoring",         lzParams.scoring);          lzParams.scoring          = NULL;
	free_score_set       ("lz.maskedScoring",   lzParams.maskedScoring);    lzParams.maskedScoring    = NULL;
	free_position_table  (targPositions);                                   targPositions             = NULL;
	free_segment_table   (anchors);                                         anchors                   = NULL;
	free_segment_table   (savedAnchors);                                    savedAnchors              = NULL;
	free_traceback       (traceback);                                       traceback                 = NULL;
	free_traceback_rows  ();
	free_if_valid        ("targCensus",         targCensus);                targCensus                = NULL;
	free_seeds           (lzParams.hitSeed);                                lzParams.hitSeed          = NULL;
	fclose_if_valid      (lzParams.capsuleFile);                            lzParams.capsuleFile      = NULL;
	free_if_valid        ("lz.capsuleFilename", lzParams.capsuleFilename);  lzParams.capsuleFilename  = NULL;
	close_capsule_file   (lzParams.capsule);                                lzParams.capsule          = NULL;
	fclose_if_valid      (lzParams.anchorsFile);                            lzParams.anchorsFile      = NULL;
	free_if_valid        ("lz.anchorsFilename", lzParams.anchorsFilename);  lzParams.anchorsFilename  = NULL;
	free_seeds           (lzParams.innerSeed);                              lzParams.innerSeed        = NULL;
	free_if_valid        ("lz.outputInfo",      lzParams.outputInfo);       lzParams.outputInfo       = NULL;
	fclose_if_valid      (lzParams.censusFile);                             lzParams.censusFile       = NULL;
	free_if_valid        ("lz.censusFilename",  lzParams.censusFilename);   lzParams.censusFilename   = NULL;
	fclose_if_valid      (lzParams.statsFile);                              lzParams.statsFile        = NULL;
	free_seed_hit_search ();
	free_quantum_search  ();

	free_if_valid        ("iz.inferFilename",   izParams.ic.inferFilename); izParams.ic.inferFilename = NULL;
	free_seeds           (izParams.hitSeed);                                izParams.hitSeed          = NULL;
	free_score_set       ("iz.scoring",         izParams.scoring);          izParams.scoring          = NULL;
	free_score_set       ("iz.maskedScoring",   izParams.maskedScoring);    izParams.maskedScoring    = NULL;

#endif // trackMemoryUsage or valgrindMemoryCheck

	// report timing stats

	dbg_timing_add (debugClockTotal);
	dbg_timing_add (debugClockQueryTotal);

	dbg_timing_report (debugClockTotal,         "total run time");
	dbg_timing_report (debugClockSeq1,          "sequence 1 I/O");
	dbg_timing_report (debugClockPosTable,      "seed position table");
	dbg_timing_report (debugClockSeq2,          "sequence 2 I/O");
	dbg_timing_report (debugClockSegTable,      "seed hit search");
	dbg_timing_report (debugClockChaining,      "chaining");
	dbg_timing_report (debugClockGappedExtend,  "gapped extension");
	dbg_timing_report (debugClockInterpolation, "interpolation");
	dbg_timing_report (debugClockOutput,        "output");
	dbg_timing_report (debugClockQueryTotal,    "total query time");

#ifdef dbgTiming
	{
	float perQuery;

	perQuery =  ((float) debugClockQueryTotal) / numQueries;
	perQuery /= clocksPerSec;

	fprintf (stderr, "%-26s %d\n",
	                 "queries:", numQueries);
	fprintf (stderr, "%-26s %.3f (%.1f per second)\n",
	                 "per query (with I/O):", perQuery, 1/perQuery);

	debugClockQueryTotal -= debugClockSeq2;
	perQuery =  ((float) debugClockQueryTotal) / numQueries;
	perQuery /= clocksPerSec;

	fprintf (stderr, "%-26s %.3f (%.1f per second)\n",
	                 "per query (w/o input):", perQuery, 1/perQuery);
	}
#endif // dbgTiming

#ifdef dbgTimingGappedExtend
	gapped_extend_timing_report (stderr);
#endif // dbgTimingGappedExtend

	return EXIT_SUCCESS;
	}

//----------
//
// capsule_target--
//  Hook up the target sequence from a capsule.
//
//----------
//
// Arguments:
//	capinfo*	cap:		The capsule info record.
//	u8**		targetRev:	Place to return a pointer to the reverse sequence.
//
// Returns:
//	A pointer to the sequence;  failures result in fatality.  The caller must
//	eventually de-allocate this by calling free_sequence().
//
//----------

static seq* capsule_target
   (capinfo*		cap,
	u8**			_targetRev)
	{
	seq*			target;
	u8*				fwd, *rev;
	char*			name;
	capseqinfo*		info;
	cappartition*	partitions;
	char*			namePool;
	u64				fwdSize, revSize, nameSize, infoSize, poolSize;
	u64				partExpected, partSize;
	u32				ix;

	// locate the mapped forward sequence

	fwd = locate_capsule_data (cap, cap_seqForward, NULL, &fwdSize);

	if (fwd == NULL)
		suicide ("bad capsule file (missing sequence)");

	if (fwdSize == 0)
		suicide ("bad capsule file, sequence length is zero");
	if (fwdSize != (unspos) fwdSize)
		suicidef ("bad capsule file, sequence length too large (0x%s)",
				  hex_64_string(fwdSize));

	if (fwd[fwdSize-1] != 0)
		suicidef ("bad capsule file, sequence not properly terminated (0x2X)",
				  fwd[fwdSize-1]);

	// locate the mapped reverse sequence

	rev = NULL; // (placate compiler)
	if (_targetRev != NULL)
		{
		rev = locate_capsule_data (cap, cap_seqReverse, NULL, &revSize);

		if (rev == NULL)
			suicide ("bad capsule file (missing reverse sequence)");
		if (revSize != fwdSize)
			suicidef ("bad capsule file, sequence lengths disagree (forward 0x%s, reverse 0x%s)",
			          hex_64_string(fwdSize), hex_64_string(revSize));

		if (rev[fwdSize-1] != 0)
			suicidef ("bad capsule file, reverse sequence not properly terminated (0x2X)",
					  rev[fwdSize-1]);
		}

	// locate the mapped name

	name = locate_capsule_data (cap, cap_seqName, NULL, &nameSize);

	if (name != NULL)
		{
		if (name[nameSize-1] != 0)
			suicidef ("bad capsule file, sequence name not properly terminated (0x2X)",
					  name[nameSize-1]);
		}

	// locate the mapped sequence info

	info = locate_capsule_data (cap, cap_seqInfo, NULL, &infoSize);

	if (info == NULL)
		suicide ("bad capsule file (missing sequence info)");

	if (infoSize != sizeof(capseqinfo))
		suicidef ("bad capsule file sequence info (expected size 0x%s, actual 0x%s)",
				  hex_64_string(sizeof(capseqinfo)), hex_64_string(infoSize));
	if (info->start == 0)
		suicidef ("bad capsule file sequence info (start = 0)");
	if (info->contig == 0)
		suicidef ("bad capsule file sequence info (contig number = 0)");
	if ((info->revCompFlags & (~rcf_revcomp)) != 0)
		suicidef ("bad capsule file sequence info (rev comp flags = %s)",
				  hex_64_string(sizeof(info->revCompFlags)));

	// locate the partition info, if needed

	partitions = NULL;
	namePool   = NULL;
	if (info->numPartitions != 0)
		{
		partExpected = ((u64) (info->numPartitions+1)) * sizeof(cappartition);

		partitions = locate_capsule_data (cap, cap_partitions, NULL, &partSize);

		if (partitions == NULL)
			suicide ("bad capsule file (missing sequence partitions)");

		if (partSize != partExpected)
			suicidef ("bad capsule file, paritions[] length mismatch (expected 0x%s, actual 0x%s)",
			          hex_64_string(partExpected), hex_64_string(partSize));

		// locate partition names

		namePool = locate_capsule_data (cap, cap_partitionNames, NULL, &poolSize);

		if (namePool == NULL)
			suicide ("bad capsule file (missing sequence partition names)");

		for (ix=0 ; ix<info->numPartitions ; ix++)
			{
			if (partitions[ix].header >= poolSize)
				suicidef ("bad capsule file, paritionName[%d] beyond array (0x%s >= 0x%s)",
				          ix,
				          hex_64_string(partitions[ix].header),
				          hex_64_string(poolSize));
			}
		}

	// create a new sequence record and hook up the mapped data

	target = new_sequence (seqposInfinity);

	target->v           = fwd;
	target->vOwner      = false;
	target->size        = fwdSize;
	target->len         = fwdSize-1;
	target->header      = name;
	target->shortHeader = name;

	if (info != NULL)
		{
		target->start        = info->start;
		target->trueLen      = info->trueLen;
		target->revCompFlags = info->revCompFlags;
		target->contig       = info->contig;
		}
	else
		{
		target->start        = 1;
		target->trueLen      = fwdSize-1;
		target->revCompFlags = rcf_forward;
		target->contig       = 1;
		}

	// hook up the partition info, if needed

	if (info->numPartitions != 0)
		{
		seqpartition* sp = &target->partition;

		sp->p         = (partition*) partitions;
		sp->size      = info->numPartitions + 1;
		sp->len       = info->numPartitions;
		sp->pool      = namePool;
		sp->poolSize  = poolSize;
		sp->poolLen   = poolSize;
		sp->poolOwner = false;
		sp->state     = seqpart_ready;
		}

	// success!

	if (_targetRev != NULL) *_targetRev = rev;
	return target;
	}

//----------
//
// capsule_position_table--
//  Hook up the target seed word position table from a capsule.
//
//----------
//
// Arguments:
//	capinfo*	cap:		The capsule info record.
//	seq*		seq:		The sequence the position table in built for.
//	seed*		hitSeed:	The seed-word the table is based on.
//	u32			step:		The step size the table is based on.
//
// Returns:
//	A pointer to the position table;  failures result in fatality.  The caller
//	must eventually de-allocate this by calling free_position_table().
//
//----------

static postable* capsule_position_table
   (capinfo*	cap,
	seq*		seq,
	seed*		hitSeed,
	u32			step)
	{
	postable*	pt;
	unspos		prevEntries;
	unspos*		last, *prev;
	u32*		asBits;
	u64			lastExpected, prevExpected, bitsExpected;
	u64			lastSize,     prevSize,     bitsSize;

	if (sizeof(unspos) != sizeof(u32))
		suicide ("internal error, capsule expects positions to be 32 bits");

	// figure out how many bytes to expect

	lastExpected = (((u64) 1) << hitSeed->weight) * sizeof(unspos);
	prevEntries  = 1 + (seq->len / step);
	prevExpected = ((u64) prevEntries) * sizeof(unspos);

	// locate the mapped last[] array

	last = locate_capsule_data (cap, cap_lastPosTable, NULL, &lastSize);

	if (last == NULL)
		suicide ("bad capsule file (missing last[] array)");

	if (lastSize != lastExpected)
		suicidef ("bad capsule file, last[] length mismatch (expected 0x%s, actual 0x%s)",
		          hex_64_string(lastExpected), hex_64_string(lastSize));

	// locate the mapped prev[] array

	prev = locate_capsule_data (cap, cap_prevPosTable, NULL, &prevSize);

	if (prev == NULL)
		suicide ("bad capsule file (missing prev[] array)");

	if (prevSize != prevExpected)
		suicidef ("bad capsule file, prev[] length mismatch (expected 0x%s, actual 0x%s)",
		          hex_64_string(prevExpected), hex_64_string(prevSize));

	// locate the mapped sequence bits

	asBits = NULL;
	if (hitSeed->type == 'R')
		{
		asBits = locate_capsule_data (cap, cap_seqBits, NULL, &bitsSize);

		if (asBits == NULL)
			suicide ("bad capsule file (missing sequence bits[] array)");

		bitsExpected = round_up_16((seq->len+3) / 4);
		if (bitsSize != bitsExpected)
			suicidef ("bad capsule file, sequence bits[] length mismatch (expected 0x%s, actual 0x%s)",
			          hex_64_string(bitsExpected), hex_64_string(bitsSize));
		}

	// create a new position table record and hook up the mapped data

	pt = new_position_table (hitSeed->weight, 0, seq->len, step,
	                         false, false, false);

	pt->last   = last;
	pt->prev   = prev;
	pt->asBits = asBits;

	// success!

	return pt;
	}

//----------
//
// set_up_hit_processor--
//  Set up variables that select and control the appropriate seed hit processor
//	function.
//
//----------
//
// Arguments:
//	control*		params:			Parameter set controlling the desired
//									.. search.
//	int				collectingCensus: true => caller will be collecting a
//									          .. census.
//	hitprocessor*	hitProc:		Place to return the function to call for
//									.. each hit to determine if it is 'good
//									.. enough'.
//	void**			hitProcInfo:	Place to return a value (usually a pointer
//									.. to some data) to pass thru with each
//									.. call to hitProc.
//
// Returns:
//  (nothing)
//
//----------
//
// Notes:
//	This routine is NOT reentrant, since some of the control data returned in
//	(*hitProcInfo) is stored in static variables that are private to this
//	routine.
//
//----------

void set_up_hit_processor
   (control*		params,
	int				collectingCensus,
	hitprocessor*	hitProc,
	void**			hitProcInfo)
	{
	static hitprocsimple	simpleInfo;
	static hitproctwin		twinInfo;
	static hitrepgappily	gappilyInfo;
	hitprocinfo*			hpInfo;
	int						filtering;

	// decide which hit processor to use when we discover a seed hit

	if (params->twinMinSpan <= 0)
		{
		if ((params->gfExtend == gfexNoExtend) && (!params->gappedExtend))
			(*hitProc) = process_for_plain_hit;
		else if (params->basicHitType == hitRecover)
			(*hitProc) = process_for_recoverable_hit;
		else
			(*hitProc) = process_for_simple_hit;
		(*hitProcInfo) = (void*) &simpleInfo;
		hpInfo         = &simpleInfo.hp;
		}
	else
		{
		(*hitProc)       = process_for_twin_hit;
		(*hitProcInfo)   = (void*) &twinInfo;
		hpInfo           = &twinInfo.hp;
		twinInfo.minSpan = params->twinMinSpan;
		twinInfo.maxSpan = params->twinMaxSpan;
		}

	params->mergeAnchors = (((*hitProc) != process_for_plain_hit)
						 && ((*hitProc) != process_for_simple_hit));

	filtering =  ((currParams->minIdentity   > 0)
	          ||  (currParams->maxIdentity   < 1)
	          ||  (currParams->minCoverage   > 0)
	          ||  (currParams->maxCoverage   < 1)
	          ||  (currParams->minContinuity > 0)
	          ||  (currParams->maxContinuity < 1)
	          ||  (currParams->minMatchCount > 0));

	// decide how to control that hit processor

	hpInfo->reporter     = collect_hsps;
	hpInfo->reporterInfo = NULL;

	if ((anchors == NULL)
	 || ((params->hspThreshold.t =='S')
	  && (!params->gappedExtend)
	  && (!params->chain)
	  && (!params->mergeAnchors)
#ifdef densityFiltering
	  && (params->maxDensity == 0)
#endif // densityFiltering
	  && (!dbgSortAnchorsByDiag)
	  && (!collectingCensus)
	  && (!filtering)))
		hpInfo->reporter = report_hsps;

	if (currParams->hspToGapped)
		{
		hpInfo->reporter           = gappily_extend_hsps;
		hpInfo->reporterInfo       = &gappilyInfo;

		gappilyInfo.seq1           = params->seq1;
		gappilyInfo.rev1           = params->rev1;
		gappilyInfo.seq2           = NULL; // (can't set this yet)
		gappilyInfo.rev2           = NULL;
		gappilyInfo.scoring        = params->scoring;
		gappilyInfo.yDrop          = params->yDrop;
		gappilyInfo.trimToPeak     = (params->yDropUntrimmed == false);
		gappilyInfo.scoreThresh    = params->gappedThreshold;
		gappilyInfo.traceback      = params->traceback;
		gappilyInfo.minIdentity    = params->minIdentity;
		gappilyInfo.maxIdentity    = params->maxIdentity;
		gappilyInfo.minCoverage    = params->minCoverage;
		gappilyInfo.maxCoverage    = params->maxCoverage;
		gappilyInfo.minContinuity  = params->minContinuity;
		gappilyInfo.maxContinuity  = params->maxContinuity;
		gappilyInfo.minMatchCount  = params->minMatchCount;
		gappilyInfo.deGapifyOutput = params->deGapifyOutput;
		}

	if (params->minMatches >= 0)
		{
		hpInfo->minMatches       = params->minMatches;
		hpInfo->maxTransversions = params->maxTransversions;
		hpInfo->filterPattern    = NULL;
		hpInfo->charToBits       = params->charToBits;
		if (params->filterCaresOnly)
			hpInfo->filterPattern = params->hitSeed->pattern;
		}
	else
		{
		hpInfo->minMatches       = -1; // (no filtering)
		hpInfo->charToBits       = params->charToBits;
		}

	if (params->gfExtend == gfexNoExtend)
		{
		hpInfo->gfExtend         = gfexNoExtend;
		hpInfo->seq1             = params->seq1;	// (sequences may be
		hpInfo->seq2             = params->seq2;	//  .. needed for filtering)
		}
	else if ((params->gfExtend == gfexExact)
	      || ((params->gfExtend >= gfexMismatch_min)
	       && (params->gfExtend <= gfexMismatch_max)))
		{
		hpInfo->gfExtend         = params->gfExtend;
		hpInfo->seq1             = params->seq1;
		hpInfo->seq2             = params->seq2;
		hpInfo->hspThreshold     = params->hspThreshold;
		hpInfo->anchors          = &anchors;
		seed_search_set_stat(isHspSearch,true);
		}
	else // if (params->gfExtend == gfexXDrop)
		{
		hpInfo->gfExtend         = gfexXDrop;
		hpInfo->seq1             = params->seq1;
		hpInfo->seq2             = params->seq2;
		hpInfo->scoring          = params->maskedScoring;
		hpInfo->xDrop            = params->xDrop;
		hpInfo->hspThreshold     = params->hspThreshold;
		hpInfo->hspZeroThreshold = (params->hspThreshold.t !='S')? 0
		                         : (params->hspThreshold.s >  0 )? params->hspThreshold.s
		                                                         : 0;
		hpInfo->anchors          = &anchors;
		hpInfo->entropicHsp      = params->entropicHsp;
		hpInfo->reportEntropy    = params->reportEntropy;
		seed_search_set_stat(isHspSearch,true);

		if (infer_scores_dbgShowIdentity)
			{
			printf ("hit_processor xDrop = "      scoreFmtSimple "\n",
			        hpInfo->xDrop);
			printf ("hit_processor hspThreshold = %s\n",
			        score_thresh_to_string (&hpInfo->hspThreshold));
			}
		}
	}

//----------
//
// start_one_strand--
//  Start alignment upon one query strand.
//
//----------
//
// Arguments:
//	seq*			target:			The sequence being searched.
//	postable*		targPositions:	A table of positions of words in target.
//	seq*			query:			The sequence(s) being searched for.
//	int				emptyAnchors:	true => clear the anchors table before
//									        .. starting.
//	hitprocessor	hitProc:		Function to call for each hit to determine
//									.. if it is 'good enough'.
//	void*			hitProcInfo:	A value (usually a pointer to some data) to
//									.. pass thru with each call to hitProc.
//
// Returns:
//  (nothing)
//
//----------

void start_one_strand
   (seq*			target,
	postable*		targPositions,
	seq*			query,
	int				emptyAnchors,
	hitprocessor	hitProc,
	void*			hitProcInfo)
	{
	u64				basesHit;

	init_output_for_strand ();

	// if we're to read achors from a file, do so

	if (currParams->anchorsFile != NULL)
		{
		if (emptyAnchors)
			empty_segment_table (anchors);
		anchors = read_segment_table
                     (currParams->anchorsFile, currParams->anchorsFilename,
		              anchors, target, query);
		return;
		}

	// find seed hits;  depending on hitProc, we will may extend these to HSPs;
	// and depending on format parameters and hitReporter, we will either report
	// these directly to the output, or we will collect them in anchors[]

	dbg_timing_sub (debugClockSegTable);

	if ((emptyAnchors) && (anchors != NULL))
		empty_segment_table (anchors);

	if (query->fileType == seq_type_qdna)
		quantum_seed_hit_search (target, targPositions,
		                         query, 0, query->len,
		                         currParams->upperCharToBits, currParams->hitSeed,
		                         currParams->maskedScoring, currParams->ballScore,
		                         hitProc, hitProcInfo);
	else
		{
		basesHit = seed_hit_search (target, targPositions,
		                            query, 0, query->len, currParams->selfCompare,
		                            currParams->upperCharToBits, currParams->hitSeed,
		                            currParams->searchLimit,
#ifdef densityFiltering
		                            currParams->maxDensity,
#endif // densityFiltering
		                            hitProc, hitProcInfo);
#ifdef densityFiltering
		if (basesHit == u64max) // maxDensity has been exceeded
			empty_segment_table (anchors);
#endif // densityFiltering
		}

	dbg_timing_add (debugClockSegTable);
	}

//----------
//
// finish_one_strand--
//  Finish alignment upon one query strand.
//
//----------
//
// Arguments:
//	seq*			target:			The sequence being searched.
//	postable*		targPositions:	A table of positions of words in target.
//	u8*				targetRev:		The reverse (NOT reverse complement) of the
//									.. target sequence, as a zero-terminated
//									.. string;  this may be NULL if the caller
//									.. doesn't need/want to supply it.  It is
//									.. only needed if we will be doing a gapped
//									.. extension.
//	seq*			query:			The sequence(s) being searched for.
//	u8*				queryRev:		The reverse (NOT reverse complement) of the
//									.. query sequence (analagous to targetRev)
//	tback*			traceback:		Memory in which to track gapped alignment
//									.. traceback.
//	census*			targCensus:		Census array for target sequence.  If this
//									.. is non-NULL, we count how many times each
//									.. target base is aligned.  This information
//									.. is used to mask positions in the target
//									.. sequence if currParams->dynamicMasking > 0.
//
// Returns:
//  (nothing)
//
//----------

//иии need to check if the use of the scoring matrix in here is copacetic with inference

void finish_one_strand
   (seq*			target,
	u8*				_targetRev,
	postable*		targPositions,
	seq*			query,
	u8*				_queryRev,
	tback*			traceback,
	census*			targCensus)
	{
	u8*				targetRev = _targetRev;
	u8*				queryRev  = _queryRev;
	alignel*		alignList = NULL;
	int				hspsAreAdaptive;
	score			lowAnchorScore = 0;

	hspsAreAdaptive = (currParams->hspThreshold.t != 'S');
	if (hspsAreAdaptive)
		{
		lowAnchorScore = anchors->lowScore;
		if ((savedAnchors != NULL)
		 && (savedAnchors->lowScore < lowAnchorScore))
			lowAnchorScore = savedAnchors->lowScore;
		}

	if ((anchors != NULL)				// merging may be necessary because
	 && (currParams->mergeAnchors))		// .. the diag hash technique used in
		{								// .. the seed search may result in
		merge_segments (anchors);		// .. duplicate or overlapping HSPs
		//fprintf (stderr, "segments for %s %c\n",
		//                   (query->partition.p != NULL)? "(partitioned query)"
		//                 : (query->useFullNames)       ? query->header
		//                                               : query->shortHeader,
		//                 ((query->revCompFlags & rcf_rev) != 0)? '-' : '+');
		//write_segments (stderr, anchors, target, query);
		}

	if (dbgSortAnchorsByDiag)
		sort_segments (anchors, qSegmentsByDiag);

	// filter HSPs by identity and/or coverage

	if ((anchors != NULL)
	 && (!currParams->gappedExtend))
		{
		if ((currParams->minIdentity > 0) || (currParams->maxIdentity < 1))
			filter_segments_by_identity (target, query, anchors,
			                             currParams->minIdentity,
			                             currParams->maxIdentity);

		if ((currParams->minCoverage > 0) || (currParams->maxCoverage < 1))
			filter_segments_by_coverage (target, query, anchors,
			                             currParams->minCoverage,
			                             currParams->maxCoverage);

		if (currParams->minMatchCount > 0)
			filter_segments_by_match_count (target, query, anchors,
			                                currParams->minMatchCount);
		}

	// if we have scoreless anchors, and we need scores, score 'em

	if ((!anchors->haveScores)
	 && ((currParams->chain) || (currParams->gappedExtend)))
		score_segments (anchors, target, query, currParams->maskedScoring);

	// reduce the set of HSPs to the best syntenic subset

	if (currParams->chain)
		{
		u32 originalNumAnchors = anchors->len;

		dbg_timing_sub  (debugClockChaining);
		reduce_to_chain (anchors, currParams->chainDiag, currParams->chainAnti,
		                 chainScale, chain_connect_penalty);
		sort_segments   (anchors, qSegmentsByPos1);
		dbg_timing_add  (debugClockChaining);

		if (dbgShowAnchors)
			fprintf (stderr, "(chaining reduced %u anchors to %u)\n",
			                   originalNumAnchors, anchors->len);
		}

	// report the set of HSPs if we don't plan to do gapped extension

	if ((anchors != NULL) && (!currParams->gappedExtend))
		{
		u32			ix;
		segment*	seg;

		for (ix=0,seg=anchors->seg ; ix<anchors->len ; ix++,seg++)
			print_match (seg->pos1, seg->pos2, seg->length, seg->s);
		}

	// if we don't plan to do gapped extension, perform dynamic masking;  note
	// that if currParams->dynamicMasking == 0, we don't actually mask, we just
	// count for the census;  further note that we don't mask the reverse
	// sequence unless the caller provided it

	if ((targCensus != NULL) && (anchors != NULL) && (!currParams->gappedExtend))
		{
		unspos numMasked;
		numMasked = census_mask_segments
		               (anchors, target->v, _targetRev, targCensus,
		                remove_interval_seeds, targPositions);
		print_x_stanza (numMasked);
		if (dbgMasking) print_m_stanza (targCensus);
		}

	// extend the HSPs to gapped alignments

	if (currParams->gappedExtend)
		{
		sthresh gappedThreshold;

		dbg_timing_sub (debugClockGappedExtend);
		if (targetRev == NULL)
			targetRev = (u8*) copy_reverse_of_string ((char*) target->v, target->len);
		if (queryRev == NULL)
			queryRev =  (u8*) copy_reverse_of_string ((char*) query->v,  query->len);

		reduce_to_points (target, query, currParams->scoring, anchors);
		if ((anchors != NULL) && (dbgShowAnchors))
			write_segments (stderr, anchors, target, query);

		gappedThreshold = currParams->gappedThreshold;
		if ((gappedThreshold.t != 'S') && (hspsAreAdaptive))
			{
			gappedThreshold.t = 'S';
			gappedThreshold.s = lowAnchorScore;
			//fprintf (stderr, "gapped threshold <- " scoreFmtSimple "\n", lowAnchorScore);
			}
		alignList = gapped_extend (target, targetRev, query, queryRev,
		                           currParams->inhibitTrivial,
		                           currParams->scoring, anchors, traceback,
		                           currParams->gappedAllBounds,
		                           currParams->yDrop,
		                           (currParams->yDropUntrimmed == false),
		                           gappedThreshold);
		dbg_timing_add (debugClockGappedExtend);
		}

	// filter gapped alignments by identity, coverage, and/or continuity

	if (alignList != NULL)
		{
		if ((currParams->minIdentity > 0) || (currParams->maxIdentity < 1))
			alignList = filter_aligns_by_identity (target, query, alignList,
			                                       currParams->minIdentity,
			                                       currParams->maxIdentity);
		if ((currParams->minCoverage > 0) || (currParams->maxCoverage < 1))
			alignList = filter_aligns_by_coverage (target, query, alignList,
			                                       currParams->minCoverage,
			                                       currParams->maxCoverage);
		if ((currParams->minContinuity > 0) || (currParams->maxContinuity < 1))
			alignList = filter_aligns_by_continuity (alignList,
			                                         currParams->minContinuity,
			                                         currParams->maxContinuity);
		if (currParams->minMatchCount > 0)
			alignList = filter_aligns_by_match_count (target, query, alignList,
			                                          currParams->minMatchCount);
		}

	// interpolate between the gapped alignments

	if ((alignList != NULL) && (currParams->innerThreshold > 0))
		{
		dbg_timing_sub (debugClockInterpolation);
		alignList = tweener_interpolate
		              (alignList, target, query,
		               currParams->selfCompare,currParams->inhibitTrivial,
		               currParams->upperCharToBits, currParams->innerSeed,
		               currParams->scoring, currParams->maskedScoring, traceback,
		               currParams->xDrop, currParams->gappedAllBounds,
		               currParams->yDrop, (currParams->yDropUntrimmed == false),
		               currParams->innerThreshold,
				       currParams->chainDiag, currParams->chainAnti,
				       chainScale, chain_connect_penalty, currParams->innerWindow);
		dbg_timing_add (debugClockInterpolation);
		}

	// print the gapped alignments

	if (alignList != NULL)
		{
		dbg_timing_sub (debugClockOutput);
		if (currParams->deGapifyOutput) print_align_list_segments (alignList);
		                           else print_align_list          (alignList);
		fflush (currParams->outputFile);
		dbg_timing_add (debugClockOutput);
		}

	// perform dynamic masking;  note that if currParams->dynamicMasking == 0, we
	// don't actually mask, we just count for the census;  further note that
	// we don't mask the reverse sequence unless the caller provided it

	if ((targCensus != NULL) && (alignList != NULL))
		{
		unspos numMasked;
		numMasked = census_mask_aligns
		               (alignList, target->v, _targetRev, targCensus,
		                remove_interval_seeds, targPositions);
		print_x_stanza (numMasked);
		if (dbgMasking) print_m_stanza (targCensus);
		}

	// cleanup

	if (alignList != NULL)
		free_align_list  (alignList);

	if ((_targetRev == NULL) && (targetRev != NULL))
		free_if_valid ("process_one_strand (targetRev)", targetRev);
	if ((_queryRev  == NULL) && (queryRev  != NULL))
		free_if_valid ("process_one_strand (queryRev)", queryRev);
	}

//----------
//
// split_anchors, swap_anchors--
//	Select single-strand anchors from a list containing anhors for both strands.
//
// Under certain configurations (e.g. for adaptive-K), collect_hsps collects
// HSPs from both strands into a single table (anchors).  Split_anchors removes
// the forward stand HSPs and moves them to a second table (savedAnchors).
// Swap_anchors switches the tables so that the forward strand HSPs are in the
// main table (anchors).
//
//----------
//
// Arguments:
//	int		id:		The id of segments on the reverse strand.
//	(none)
//
// Returns:
//  (nothing)
//
//----------

void split_anchors (int id)
	{
	if (savedAnchors == NULL)
		savedAnchors = new_segment_table (numDefaultAnchors,
		                                  anchors->coverageLimit);
	else
		empty_segment_table (savedAnchors);

	split_segment_table (anchors, id, &savedAnchors);

	//printf ("\nanchors:\n");
	//dump_segments (stdout, anchors, NULL, NULL);
	//printf ("\nleftovers:\n");
	//dump_segments (stdout, savedAnchors, NULL, NULL);
	//printf ("\n");
	}


void swap_anchors (void)
	{ segtable* a = savedAnchors;  savedAnchors = anchors;  anchors = a; }

//----------
// [[-- a chain connection penalty function --]]
//
// chain_connect_penalty--
//	Compute penalty for connecting two segments in the chain.
//
// Arguments and Return value: (see chain.h)
//
// Note bene:
//	x = pos1, y = pos2, diag = x-y
//	diag increases toward southeast
//
//----------

#define debugChaining_1                                                      \
	if (chain_dbgChaining)                                                   \
		fprintf (stderr,                                                     \
		         "    diagDiff=" sgnposFmt " numSubs=" sgnposFmt "\n",       \
				 diagDiff, numSubs);

#define debugChaining_2                                                      \
	if (chain_dbgChaining)                                                   \
		fprintf (stderr, "    base_penalty=%.2f\n", penalty);

#define debugChaining_3                                                      \
	if (chain_dbgChaining)                                                   \
		fprintf (stderr, "    penalty=%.2f\n", penalty);

#define debugChaining_4                                                      \
	if (chain_dbgChaining)                                                   \
		fprintf (stderr, "    penalty=%.2f (" scoreFmtSimple ")\n",          \
		                 penalty, currParams->scoring->sub[rCh][cCh]);

#define debugChaining_5                                                      \
	if (chain_dbgChaining)                                                   \
		{                                                                    \
		if (penalty > bestPossibleScore)                                     \
			fprintf (stderr, "    returning " scoreFmtSimple "\n",           \
			                 bestPossibleScore);                             \
		else                                                                 \
			fprintf (stderr, "    returning " scoreFmtSimple "\n",           \
			                 (score) penalty);                               \
		}


static score chain_connect_penalty
   (segment*	seg1,
	segment*	seg2,
	int			scale)
	{
	unspos		xEnd, yEnd;
	sgnpos		diag1, diag2, diagDiff;
	sgnpos		numSubs;	// number of substitutions needed to get from end
	double		penalty;	// .. of segment 1 to beginning of segment 2
	u8			rCh, cCh;

	if ((seg2->pos1 <= seg1->pos1) || (seg2->pos2 <= seg1->pos2))
		suicide ("HSPs improperly ordered for chaining");

	xEnd  = seg1->pos1 + seg1->length - 1;
	yEnd  = seg1->pos2 + seg1->length - 1;

	diag1 = diagNumber (seg1->pos1, seg1->pos2);
	diag2 = diagNumber (seg2->pos1, seg2->pos2);

	diagDiff = diag2 - diag1;
	if (diagDiff >= 0)
		{					// segment 1's diagonal is above segment 2's
		numSubs = ((sgnpos) seg2->pos2) - ((sgnpos) yEnd) - 1;
		}
	else
		{					// segment 1's diagonal is below segment 2's
		numSubs  = ((sgnpos) seg2->pos1) - ((sgnpos) xEnd) - 1;
		diagDiff = -diagDiff;
		}

	debugChaining_1;

	// nota bene: penalty is declared as double to allow it to overflow the
	//            regular score type;  after we compute the penalty we clip it
	//            to worst penalty as we return it

	penalty = diagDiff * currParams->chainDiag;
	debugChaining_2;
	if (numSubs >= 0)
		{
		penalty +=   numSubs  * currParams->chainAnti;
		debugChaining_3;
		}
	else
		{
		rCh = currParams->scoring->rowChars[0];
		cCh = currParams->scoring->colChars[0];
		penalty += (-numSubs) * scale * currParams->scoring->sub[rCh][cCh];
		debugChaining_4;
		}

	debugChaining_5;
	if (penalty > bestPossibleScore) return bestPossibleScore;
	                            else return (score) penalty;
	}

//----------
// [[-- a report_census_intervals callback function --]]
//
// remove_interval_seeds--
//	Remove seeds from the target sequence position table that are about to be
//	rendered meaningless by dynamic masking.
//
//----------
//
// Arguments:
//	unspos b, e:	The interval, in the target sequence, that is about to be
//					.. masked.  Origin-1, inclusive.
//	void* info:		(really postable*) The position table (targPositions).
//
// Returns:
//  (nothing)
//
//----------
//
// Notes:
//	(1)	Intervals given to this routine are origin-1 inclusive, while the
//		intervals it passes along to mask_seed_position_table are origin-0
//		end-exclusive.  Further, we have to expand the interval on both ends
//		by the length of the seed (minus 1).
//
//		For example, suppose the input interval is 20..40 (11 bp) and the seed
//		length is 10.  In the diagram below, * is a base in the input interval,
//		x is a base in the expanded interval, and o is any other base.  Numbers
//		on the top are origin-1; numbers on the bottom are origin-0.  The output
//		interval is 10..49.
//
//		    1        11       20                  40       49
//			v         v        v                   v        v
//			ooooooooooxxxxxxxxx*********************xxxxxxxxxoooooooo ...
//			^         ^        ^                   ^        ^
//			0        10       19                  39       48
//
//		The reason for expansion is that we want to eliminate any seed that
//		contains any * in the diagram.  The * can occur at any position in the
//		10 bp seed, so we must expand by 9 bp on each end.
//
//----------

static void remove_interval_seeds (unspos b, unspos e, void* info)
	{
	postable*	pt = (postable*) info;
	seq*		target  = currParams->seq1;
	seed*		hitSeed = currParams->hitSeed;
	u32			seedLen = (unsigned) hitSeed->length;
	const s8*	upperCharToBits = currParams->upperCharToBits;

	// adjust the interval endpoints to account for the seed length

	if (b < seedLen) b =  1;
	            else b -= seedLen - 1;

	if (e >= target->len - (seedLen-1)) e =  target->len;
	                               else e += seedLen - 1;

	// remove masked seeds from the table

	mask_seed_position_table (pt, target, b-1, e, upperCharToBits, hitSeed);
	}

//----------
// [[-- a seed hit reporter function --]]
//
// report_hsps--
//	Report a seed hit or HSP (i.e. just write it to output).
//
// Arguments and Return value: (see seed_search.h)
//
//----------

static u32 report_hsps
   (arg_dont_complain(void* info),
	unspos	pos1,
	unspos	pos2,
	unspos	length,
	score	s)
	{
	unspos	s1, s2;

	// report this hit/HSP

	print_match (pos1-length, pos2-length, length, s);

	if (dbgShowHsps)
		{
		fprintf (stderr, "\n");
		dump_aligned_nucleotides (stderr,
		                          currParams->seq1, pos1-length,
		                          currParams->seq2, pos2-length,
		                          length);
		}

	if (!currParams->mirror)
		return length;

	// report the mirror of this hit/HSP

	if (currParams->seq1->revCompFlags == currParams->seq2->revCompFlags)
		{
		s1 = pos1;
		s2 = pos2;
		}
	else
		{
		s1 = (currParams->seq1->len) - pos1 + length;
		s2 = (currParams->seq2->len) - pos2 + length;
		if ((s2 == pos1) && (s1 == pos2)) return length;
		}

	print_match (s2-length, s1-length, length, s);

	if (dbgShowHsps)
		{
		fprintf (stderr, "\n");
		dump_aligned_nucleotides (stderr,
		                          currParams->seq1, s2-length,
		                          currParams->seq2, s1-length,
		                          length);
		}

	return length;
	}

//----------
// [[-- a seed hit reporter function --]]
//
// collect_hsps--
//	Collect a seed hit or HSP.
//
// Arguments and Return value: (see seed_search.h)
//
//----------

static u32 collect_hsps
   (arg_dont_complain(void* info),
	unspos	pos1,
	unspos	pos2,
	unspos	length,
	score	s)
	{
	unspos	s1, s2;

	// add this hit/HSP to the list of anchors;  note that we use the strand
	// as the id field, so that if we happen to be collecting segments from
	// both strands, we can separate them later

	if (dbgShowAnchors)
		fprintf (stderr, "adding segment " unsposSlashFmt " " unsposFmt
		                 " diag=" sgnposFmt "\n",
		                 pos1-length, pos2-length, length,
		                 diagNumber(pos1-length,pos2-length));

	anchors = add_segment (anchors, pos1-length, pos2-length, length, s,
	                       /* id */ currParams->seq2->revCompFlags);

	if (dbgShowAnchors)
		fprintf (stderr, "(now have %u anchors)\n", anchors->len);

	if (dbgShowHsps)
		{
		fprintf (stderr, "\n");
		dump_aligned_nucleotides
			 (stderr, currParams->seq1, pos1-length, currParams->seq2, pos2-length, length);
		}

	if (!currParams->mirror)
		return length;

	// add the mirror of this hit/HSP to the list of anchors

	if (currParams->seq1->revCompFlags == currParams->seq2->revCompFlags)
		{
		s1 = pos1;
		s2 = pos2;
		}
	else
		{
		s1 = (currParams->seq1->len) - pos1 + length;
		s2 = (currParams->seq2->len) - pos2 + length;
		if ((s2 == pos1) && (s1 == pos2)) return length;
		}

	if (dbgShowAnchors)
		fprintf (stderr, "adding segment " unsposSlashFmt " " unsposFmt
		                 " diag=" sgnposFmt "\n",
		                 s2-length, s1-length, length,
		                 diagNumber(s2-length,s1-length));

	anchors = add_segment (anchors, s2-length, s1-length, length, s,
	                       /* id */ currParams->seq2->revCompFlags);

	if (dbgShowAnchors)
		fprintf (stderr, "(now have %u anchors)\n", anchors->len);

	if (dbgShowHsps)
		{
		fprintf (stderr, "\n");
		dump_aligned_nucleotides
			 (stderr, currParams->seq1, s2-length, currParams->seq2, s1-length, length);
		}

	return 2*length;
	}

//----------
//
// parse_options--
//  Parse command line options.
//
//----------
//
// Arguments:
//	argc, argv:	(as per main)
//	control*	lzParams:	Control data to fill in for the primary alignment.
//	control*	izParams:	Control data to fill in for inference alignments.
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

// pre-canned expansion arguments

typedef struct exparg
	{
	char*	argName;
	char*	expansion;
	} exparg;

exparg expanders[] =
	{
// old expansions:
//	{ "--yasra98",      "T=2 Z=20 --match=1,6 O=6 E=1 Y=20 K=22 L=30 --identity=98..100" },
//	{ "--yasra95",      "T=2 Z=20 --match=1,5 O=8 E=1 Y=40 K=40 L=40 --identity=95..100" },
//	{ "--yasra90",      "T=2 Z=20 --match=1,5 O=6 E=1 Y=40 K=40 L=40 --identity=90..100" },
//	{ "--yasra85",      "T=2      --match=1,2 O=4 E=1 Y=20 K=22 L=30 --identity=85..100" },
//	{ "--yasra75",      "T=2      --match=1,1 O=3 E=1 Y=20 K=22 L=30 --identity=75..100" }
//	{ "--yasra95short", "T=2 --match=1,7 K=10 --nogapped --noentropy --identity=95..100" },
//	{ "--yasra85short", "T=2 --match=1,2 K=11 --nogapped --noentropy --identity=85..100" }
	{ "--yasra98",      "T=2 Z=20 --match=1,6 O=8 E=1 Y=20 K=22 L=30 --identity=98..100" },
	{ "--yasra95",      "T=2 Z=20 --match=1,5 O=8 E=1 Y=20 K=22 L=30 --identity=95..100" },
	{ "--yasra90",      "T=2 Z=20 --match=1,5 O=6 E=1 Y=20 K=22 L=30 --identity=90..100" },
	{ "--yasra85",      "T=2      --match=1,2 O=4 E=1 Y=20 K=22 L=30 --identity=85..100" },
	{ "--yasra75",      "T=2      --match=1,1 O=3 E=1 Y=20 K=22 L=30 --identity=75..100" },
	{ "--yasra95short", "T=2      --match=1,7 O=6 E=1 Y=14 K=10 L=14 --identity=95..100" },
	{ "--yasra85short", "T=2      --match=1,3 O=4 E=1 Y=14 K=11 L=14 --identity=85..100" }
	};
#define numExpanders ((int)(sizeof(expanders)/sizeof(exparg)))


static void usage (void);
static void usage (void)
	{
	fprintf (helpout, "%s-- Local Alignment Search Tool, blastZ-like\n",
	                  programName);
	fprintf (helpout, "  (version %s.%s.%s released %s",
	                  programVersionMajor, programVersionMinor, programVersionSubMinor, programRevisionDate);
	if      (scoreType == 'F') fprintf (helpout, ", floating point scores");
	else if (scoreType == 'D') fprintf (helpout, ", double floating point scores");
	fprintf (helpout, ")\n");
	fprintf (helpout, "usage: %s target [query] [options]\n", programName);

	fprintf (helpout, "  (common options;  use --help for a more extensive list)\n");
	fprintf (helpout, "  target, query          specifiers or files, containing sequences to align\n");
	fprintf (helpout, "                         (use --help=files for more details)\n");

	fprintf (helpout, "  --seed=<pattern>       set seed pattern (12of19, 14of22, or general pattern)\n");
	fprintf (helpout, "                         (default is %s)\n",
	         defaultSeedString);
	fprintf (helpout, "  --[no]transition       allow (or don't) one transition in a seed hit\n");
	fprintf (helpout, "                         (by default %s)\n",
	         (defaultParams.withTrans == 0)? "the seed must match as is" :
	         (defaultParams.withTrans == 1)? "a transition is allowed"
	                                       : "two transitions are allowed");

	fprintf (helpout, "  --[no]chain            perform chaining\n");
	fprintf (helpout, "                         (by default %s)\n",
	         (!defaultParams.chain)?  "no chaining is performed"
	                                : "chaining is performed");

	fprintf (helpout, "  --[no]gapped           perform gapped alignment (instead of gap-free)\n");
	fprintf (helpout, "                         (by default %s)\n",
	         (!defaultParams.gappedExtend)? "gapped alignment is not performed"
	                                      : "gapped alignment is performed");

	fprintf (helpout, "  --step=<length>        set step length (default is %u)\n",
	         defaultParams.step);

	fprintf (helpout, "  --strand=both          search both strands\n");
	fprintf (helpout, "  --strand=plus          search + strand only (matching strand of query spec)\n");
	fprintf (helpout, "                         (by default %s)\n",
	         (defaultParams.whichStrand == 0)? "only + strand is searched" :
	         (defaultParams.whichStrand <  0)? "only - strand is searched"
	                                         : "both strands are searched");

	fprintf (helpout, "  --scores=<file>        read substitution and gap scores from a file\n");
	fprintf (helpout, "  --xdrop=<score>        set x-drop threshold (default is 10sub[A][A])\n");
	fprintf (helpout, "  --ydrop=<score>        set y-drop threshold (default is open+300extend)\n");

	fprintf (helpout, "  --infer[=<control>]    infer scores from the sequences, then use them\n");
	fprintf (helpout, "                         all inference options are read from the control file\n");

	fprintf (helpout, "  --hspthresh=<score>    set threshold for high scoring pairs (default is %s)\n",
	         score_thresh_to_string (&defaultParams.hspThreshold));
	fprintf (helpout, "                         ungapped extensions scoring lower are discarded\n");
	fprintf (helpout, "                         <score> can also be a percentage or base count\n");

	fprintf (helpout, "  --gappedthresh=<score> set threshold for gapped alignments\n");
	fprintf (helpout, "                         gapped extensions scoring lower are discarded\n");
	fprintf (helpout, "                         <score> can also be a percentage or base count\n");
	fprintf (helpout, "                         (default is to use same value as --hspthresh)\n");

	fprintf (helpout, "  --help                 list all options\n");
	fprintf (helpout, "  --help=files           list information about file specifiers\n");
	fprintf (helpout, "  --help=shortcuts       list blastz-compatible shortcuts\n");
	fprintf (helpout, "  --help=yasra           list yasra-specific shortcuts\n");

	exit (EXIT_FAILURE);
	}


static void all_options (void);
static void all_options (void)
	{
	int ix;

	fprintf (helpout, "  target[[start..end]]   spec/file containing target sequence (fasta, nib, 2bit\n");
	fprintf (helpout, "                         or hsx);  [start..end] defines a subrange of the file\n");
	fprintf (helpout, "                         (use --help=files for more details)\n");
	fprintf (helpout, "  query[[start..end]]    spec/file containing query sequences;  if absent,\n");
	fprintf (helpout, "                         queries come from stdin (if needed)\n");
	fprintf (helpout, "  --self                 the target sequence is also the query\n");

	fprintf (helpout, "  --seed=match<length>   use a word with no gaps instead of a seed pattern\n");
	fprintf (helpout, "  --seed=half<length>    use space-free half-weight word instead of seed pattern\n");
	fprintf (helpout, "  --[no]transition[=2]   allow one or two transitions in a seed hit\n");
	fprintf (helpout, "                         (by default %s)\n",
	         (defaultParams.withTrans == 0)? "the seed must match as is" :
	         (defaultParams.withTrans == 1)? "a transition is allowed"
	                                       : "two transitions are allowed");

	fprintf (helpout, "  --word=<bits>          set max bits for word hash;  use this to trade time for\n");
	fprintf (helpout, "                         memory, eliminating thrashing for heavy seeds\n");
	fprintf (helpout, "                         (default is %d bits)\n",
	         defaultParams.maxIndexBits);

	fprintf (helpout, "  --filter=<T>,<M>       filter seed hits, requiring at least M matches and\n");
	fprintf (helpout, "                         allowing no more than T transversions\n");
	if (defaultParams.minMatches < 0)
		fprintf (helpout, "                         (default is no filtering)\n");
	else if (defaultParams.maxTransversions < 0)
		fprintf (helpout, "                         (default is to require %d matches)\n",
		         defaultParams.minMatches);
	else if (defaultParams.maxTransversions == 0)
		fprintf (helpout, "                         (default is to require %d matches, no transversions)\n",
		         defaultParams.minMatches);
	else
		fprintf (helpout, "                         (default is %d matches/%d transversions)\n",
		         defaultParams.minMatches,defaultParams.maxTransversions);

	fprintf (helpout, "  --notwins              require just one seed hit\n");
	fprintf (helpout, "  --twins=<min>..<maxgap> require two nearby seed hits on the same diagonal\n");
	if (defaultTwinsYes)
		fprintf (helpout, "                         (default is twins with %d:%d bp gap)\n",
		         defaultTwinMinGap,defaultTwinMaxGap);
	else
		fprintf (helpout, "                         (default is twins aren't required)\n");

#ifndef noSeedHitQueue
	fprintf (helpout, "  --seedqueue=<entries>  set number of entries in seed hit queue\n");
	fprintf (helpout, "                         (default is %d)\n",
	                  defaultParams.seedHitQueueSize);
#endif // not noSeedHitQueue

	fprintf (helpout, "  --segments=<file>      read anchor segments from a file, instead of\n");
	fprintf (helpout, "                         discovering anchors via seeding\n");

	fprintf (helpout, "  --norecoverseeds       don't recover hash-collision seed hits\n");
	fprintf (helpout, "  --recoverseeds         recover hash-collision seed hits\n");
	if (defaultParams.basicHitType == hitRecover)
		fprintf (helpout, "                         (default is to recover seed hits)\n");
	else
		fprintf (helpout, "                         (default is not to recover seed hits)\n");

	fprintf (helpout, "  --step=<length>        set step length (default is %u)\n",
	         defaultParams.step);

	fprintf (helpout, "  --strand=both          search both strands\n");
	fprintf (helpout, "  --strand=plus          search + strand only (matching strand of query spec)\n");
	fprintf (helpout, "  --strand=minus         search - strand only (opposite strand of query spec)\n");
	fprintf (helpout, "                         (by default %s)\n",
	         (defaultParams.whichStrand == 0)? "only + strand is searched" :
	         (defaultParams.whichStrand <  0)? "only - strand is searched"
	                                         : "both strands are searched");

	fprintf (helpout, "  --ambiguous=n          treat N as an ambiguous nucleotide\n");
	fprintf (helpout, "                         (by default N is treated as a sequence splicing\n");
	fprintf (helpout, "                          character)\n");

	fprintf (helpout, "  --ambiguous=iupac      treat any ambiguous IUPAC-IUB character as a\n");
	fprintf (helpout, "                         completely ambiguous nucleotide\n");
	fprintf (helpout, "                         (by default any sequence file with B,D,H,K,M,R,S,V,W,Y\n");
	fprintf (helpout, "                          is rejected)\n");

	fprintf (helpout, "  --[no]gfextend         perform gap-free extension of seed hits to HSPs\n");
	fprintf (helpout, "                         (by default %s)\n",
	           (defaultParams.gfExtend == gfexNoExtend)?  "no extension is performed"
	         : (defaultParams.gfExtend == gfexExact)?     "exact match extension is performed"
	                                                    : "extension is performed");

	fprintf (helpout, "  --[no]chain            perform chaining\n");
	fprintf (helpout, "  --chain=<diag,anti>    perform chaining with given penalties for diagonal and\n");
	fprintf (helpout, "                         anti-diagonal\n");
	fprintf (helpout, "                         (by default %s)\n",
	         (!defaultParams.chain)?  "no chaining is performed"
	                                : "chaining is performed");

	fprintf (helpout, "  --[no]gapped           perform gapped alignment (instead of gap-free)\n");
	fprintf (helpout, "                         (by default %s)\n",
	         (!defaultParams.gappedExtend)? "gapped alignment is not performed"
	                                      : "gapped alignment is performed");
	fprintf (helpout, "  --notrivial            do not output a trivial self-alignment block if the\n");
	fprintf (helpout, "                         target and query happen to be identical\n");

	fprintf (helpout, "  --scores=<file>        read substitution scores from a file\n");
	fprintf (helpout, "                         (default is HOXD70)\n");
	fprintf (helpout, "  --match=<R>,<P>        scores are +R/-P for match/mismatch\n");
	fprintf (helpout, "  --gap=<open,extend>    set gap open and extend penalties (default is " scoreFmtSimple "," scoreFmtSimple ")\n",
	                 HOXD70_open, HOXD70_extend);
	fprintf (helpout, "  --xdrop=<score>        set x-drop threshold (default is 10*sub[A][A])\n");
	fprintf (helpout, "  --ydrop=<score>        set y-drop threshold (default is open+300extend)\n");
	fprintf (helpout, "  --noxtrim              if x-drop extension encounters end of sequence, don't\n");
	fprintf (helpout, "                         trim back to peak score (use this for short reads)\n");
	fprintf (helpout, "  --noytrim              if y-drop extension encounters end of sequence, don't\n");
	fprintf (helpout, "                         trim back to peak score (use this for short reads)\n");

	fprintf (helpout, "  --infer=<control>      infer scores from the sequences, then use them\n");
	fprintf (helpout, "  --inferonly=<control>  infer scores but don't use them (requires --infscores)\n");
	fprintf (helpout, "                         all inference options are read from the control file\n");
	fprintf (helpout, "  --infscores[=<file>]   write inferred scores to a file\n");

	fprintf (helpout, "  --hspthresh=<score>    set threshold for high scoring pairs (default is %s)\n",
	         score_thresh_to_string (&defaultParams.hspThreshold));
	fprintf (helpout, "                         ungapped extensions scoring lower are discarded\n");
	fprintf (helpout, "                         <score> can also be a percentage or base count\n");

	fprintf (helpout, "  --exact=<length>       set threshold for exact matches\n");
	fprintf (helpout, "                         if specified, exact matches are found rather than high\n");
	fprintf (helpout, "                         scoring pairs (replaces --hspthresh)\n");
	fprintf (helpout, "  --mismatch=<N>,<length> set threshold for mismatches\n");
	fprintf (helpout, "                         if specified, N-mismatch segments are found rather\n");
	fprintf (helpout, "                         than high scoring pairs (replaces --hspthresh)\n");

	fprintf (helpout, "  --inner=<score>        set threshold for HSPs during interpolation\n");
	if (defaultParams.innerThreshold <= 0)
		fprintf (helpout, "                         (default is no interpolation)\n");
	else
		fprintf (helpout, "                         (default is " scoreFmtSimple ")\n",
		         defaultParams.innerThreshold);

	fprintf (helpout, "  --gappedthresh=<score> set threshold for gapped alignments\n");
	fprintf (helpout, "                         gapped extensions scoring lower are discarded\n");
	fprintf (helpout, "                         <score> can also be a percentage or base count\n");
	fprintf (helpout, "                         (default is to use same value as --hspthresh)\n");

	fprintf (helpout, "  --ball=<score>[%%]      set minimum score required of words 'in' a quantum ball\n");

	fprintf (helpout, "  --[no]entropy          involve entropy in filtering high scoring pairs\n");
	fprintf (helpout, "                         (default is \"%s\")\n",
	         (defaultParams.entropicHsp)? "entropy"
	                                    : "noentropy");

	fprintf (helpout, "  --nomirror             don't report mirror-image alignments when using --self\n");
	fprintf (helpout, "                         (default is to skip processing them, but recreate them\n");
	fprintf (helpout, "                         in the output)\n");

	fprintf (helpout, "  --traceback=<bytes>    space for trace-back information\n");
	fprintf (helpout, "                         (default is %s)\n",
	         unitize(defaultParams.tracebackMem));

	fprintf (helpout, "  --maxwordcount=<limit>[%%] limit seed word-repeats in target\n");
	fprintf (helpout, "                         words occurring too often are not used in seed hits\n");
	fprintf (helpout, "                         (default is no word-repeat limit)\n");

	fprintf (helpout, "  --masking=<count>      mask any position in target hit this many times\n");
	fprintf (helpout, "                         zero indicates no masking\n");
	if (defaultParams.dynamicMasking <= 0)
		fprintf (helpout, "                         (default is no masking)\n");
	else
		fprintf (helpout, "                         (default is %d)\n",
		         defaultParams.dynamicMasking);

	fprintf (helpout, "  --[no]census[=<file>]  count/report how many times each target base aligns\n");
	fprintf (helpout, "                         (default is %s)\n",
	         (defaultParams.reportCensus)? "to report census"
	                                     : "to not report census");

	fprintf (helpout, "  --identity=<min>[..<max>] filter alignments by percent identity\n");
	fprintf (helpout, "                         0<=min<=max<=100;  blocks (or HSPs) outside min..max\n");
	fprintf (helpout, "                         are discarded\n");
	fprintf (helpout, "                         (default is no identity filtering)\n");

	fprintf (helpout, "  --coverage=<min>[..<max>] filter alignments by percentage of query covered\n");
	fprintf (helpout, "                         0<=min<=max<=100;  blocks (or HSPs) outside min..max\n");
	fprintf (helpout, "                         are discarded\n");
	fprintf (helpout, "                         (default is no query coverage filtering)\n");

	fprintf (helpout, "  --continuity=<min>[..<max>] filter alignments by percent continuity\n");
	fprintf (helpout, "                         0<=min<=max<=100;  blocks (or HSPs) outside min..max\n");
	fprintf (helpout, "                         are discarded\n");
	fprintf (helpout, "                         (default is no continuity filtering)\n");

	fprintf (helpout, "  --matchcount=<min> filter alignments by match-count\n");
	fprintf (helpout, "                         0<min;  blocks (or HSPs) with fewer than min matched\n");
	fprintf (helpout, "                         bases are discarded\n");
	fprintf (helpout, "                         (default is no match-count filtering)\n");

#ifdef densityFiltering
	fprintf (helpout, "  --density=<max>        filter query sequences by alignment density\n");
	fprintf (helpout, "                         sequences with bases_aligned/bases > max\n");
	fprintf (helpout, "                         are discarded\n");
	fprintf (helpout, "                         (default is no query density filtering)\n");
#endif // densityFiltering

//	fprintf (helpout, "  --[no]laj              backward compatibility for laj\n");
//	fprintf (helpout, "                         (default is %s)\n",
//	         (defaultParams.lajCompatible)? "to be backward compatible"
//	                                      : "not to bother with backward compatibility");

	fprintf (helpout, "  --output=<file>        specify output alignment file;  otherwise alignments\n");
	fprintf (helpout, "                         are written to stdout\n");
	fprintf (helpout, "  --format=<type>        specify output format; one of lav, axt, maf, cigar,\n");
	fprintf (helpout, "                         rdotplot, text or general\n");
	fprintf (helpout, "                         (use --help=formats for more details)\n");
	fprintf (helpout, "                         (by default output is %s)\n",
	         formatNames[defaultParams.outputFormat]);
	fprintf (helpout, "  --markend              Write a comment at the end of the output file\n");

	fprintf (helpout, "  --rdotplot=<file>      create an output file suitable for plotting in R.\n");

	fprintf (helpout, "  --verbosity=<level>    set info level (0 is minimum, 10 is everything)\n");
	fprintf (helpout, "                         (default is %d)\n",
	         defaultParams.verbosity);

	fprintf (helpout, "  --[no]runtime          report runtime in the output file\n");
	fprintf (helpout, "                         (default is %s)\n",
	         (defaultParams.reportTiming)? "to report runtime"
	                                     : "to not report runtime");

	fprintf (helpout, "  --tableonly[=count]    just produce the target position table, don't\n");
	fprintf (helpout, "                         search for seeds\n");

	fprintf (helpout, "  --writesegments=<file> just produce the anchor segments table, don't\n");
	fprintf (helpout, "                         perform gapped alignment\n");

	fprintf (helpout, "  --writecapsule=<file>  write the target and seed word table to a file\n");
	fprintf (helpout, "  --targetcapsule=<file> read the target seed word table from a file\n");
	fprintf (helpout, "                         (this replaces the target specifier)\n");

#ifdef collect_stats
	fprintf (helpout, "  --[no]stats[=<file>]   show search statistics (or don't)\n");
	fprintf (helpout, "                         (default is %s)\n",
	         (defaultParams.showStats)? "to show shats"
	                                  : "to not show stats");
#endif

	fprintf (helpout, "  --progress=<n>         report processing of every nth query\n");

	fprintf (helpout, "  --version              report the program version and quit\n");
	fprintf (helpout, "  --help                 list all options\n");
	fprintf (helpout, "  --help=files           list information about file specifiers\n");
	fprintf (helpout, "  --help=formats         list information about output file formats\n");
	fprintf (helpout, "  --help=shortcuts       list blastz-compatible shortcuts\n");
	fprintf (helpout, "  --help=yasra           list yasra-specific shortcuts\n");

	for (ix=0 ; ix<numExpanders ; ix++)
		{
		if (strcmp_prefix (expanders[ix].argName, "--yasra") == 0) continue;
		fprintf (helpout, "  %-22s (%s)\n",
		                  expanders[ix].argName, expanders[ix].expansion);
		}

	exit (EXIT_FAILURE);
	}

static void file_options (void);
static void file_options (void)
	{
	fprintf (helpout, "usage: %s target [query] [options]\n", programName);
	fprintf (helpout, "\n");
	fprintf (helpout, "target is required unless replaced by the --targetcapsule option.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "query is not required;  it can be replaced by --self or in some cases (e.g.\n");
	fprintf (helpout, "with --tableonly or --writecapsule) no query sequence is needed.  If a query\n");
	fprintf (helpout, "sequence is needed and the query field is absent, the sequence is read from\n");
	fprintf (helpout, "stdin.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "The general form of both target and query specifiers is\n");
	fprintf (helpout, "    [nickname::]filename[/selectname][[actions]][-]\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "The filename field is required; all other fields are optional.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "Fields:\n");
	fprintf (helpout, "  nickname            name to use for this sequence in any output files\n");
	fprintf (helpout, "  filename            file (or path) to fasta, nib, 2bit or hsx file\n");
	fprintf (helpout, "  selectname          read only a single sequence from the file\n");
	fprintf (helpout, "                      (only valid for 2bit or hsx)\n");
	fprintf (helpout, "  actions             list of pre-processing actions;  enclosed in square\n");
	fprintf (helpout, "                      brackets and comma-separated;  see list of actions below\n");
	fprintf (helpout, "  - (minus sign)      use reverse complement of the sequence\n");
	fprintf (helpout, "                      (equivalent to the revcomp action listed below)\n");
	fprintf (helpout, "\n");

	print_file_actions (helpout);

	exit (EXIT_FAILURE);
	}

static void format_options (void);
static void format_options (void)
	{
	int ix, nameLen, needComma, lineWidth;
	char* name;

	fprintf (helpout, "Lastz Output File Formats\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "(output is written to stdout unless the --output option is used)\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "LAV\n");
	fprintf (helpout, "    LAV format is the format that blastz produced, and is the default.  It\n");
	fprintf (helpout, "    reports alignment blocks grouped by 'contig' and strand, and describes the\n");
	fprintf (helpout, "    blocks by listing the coordinates of ungapped segments.  It does not display\n");
	fprintf (helpout, "    the nucleotides.  For more deatils see the lastz readme file.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "    The option --format=lav+text adds a textual display of each alignment\n");
	fprintf (helpout, "    block, intermixed with the lav format.  Such files are unlikely to be\n");
	fprintf (helpout, "    recognized by any lav-reading program.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "AXT\n");
	fprintf (helpout, "    AXT format is a pairwise alignment format.  As of Jan/2009, a spec for AXT\n");
	fprintf (helpout, "    files can be found at\n");
	fprintf (helpout, "        genome.ucsc.edu/goldenPath/help/axt.html\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "    The option --format=axt+ displays additional statistics with each block,\n");
	fprintf (helpout, "    in the form of comments.  The exact content of these comment lines may\n");
	fprintf (helpout, "    change in future releases of lastz.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "MAF\n");
	fprintf (helpout, "    MAF format is a multiple alignment format.  As of Jan/2009, a spec for MAF\n");
	fprintf (helpout, "    files can be found at\n");
	fprintf (helpout, "        genome.ucsc.edu/FAQ/FAQformat#format5\n");
	fprintf (helpout, "    The MAF files produced by lastz have exactly two sequences per block.  The\n");
	fprintf (helpout, "    first sequence always comes from the target sequence file, the second from\n");
	fprintf (helpout, "    the query.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "    The option --format=maf+ displays additional statistics with each block,\n");
	fprintf (helpout, "    in the form of comments.  The exact content of these comment lines may\n");
	fprintf (helpout, "    change in future releases of lastz.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "    The option --format=maf- inhibits the maf header and any comments.  This\n");
	fprintf (helpout, "    makes it suitable for catenating output from multiple runs.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "SAM\n");
	fprintf (helpout, "    SAM format is a pairwise alignment format used primarily for short-read\n");
	fprintf (helpout, "    mapping.  It is imperative that the query sequence(s) be short reads.  By\n");
	fprintf (helpout, "    default \"hard clipping\" is used when alignments don't reach the end of a\n");
	fprintf (helpout, "    query (see the SAM spec for what that means).  The option --format=softsam\n");
	fprintf (helpout, "    will use \"soft clipping\" instead.  As of Oct/2009, a spec for SAM files\n");
	fprintf (helpout, "    can be found at\n");
	fprintf (helpout, "        samtools.sourceforge.net/SAM1.pdf\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "    The option --format=sam- inhibits the sam header lines.  This makes it\n");
	fprintf (helpout, "    suitable for catenating output from multiple runs.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "CIGAR\n");
	fprintf (helpout, "    CIGAR format is a pairwise alignment format that describes alignment blocks\n");
	fprintf (helpout, "    in a run-length format.  As of Jan/2009, a spec for CIGAR files can be\n");
	fprintf (helpout, "    found at\n");
	fprintf (helpout, "        may2005.archive.ensembl.org/Docs/wiki/html/EnsemblDocs/CigarFormat.html\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "segments\n");
	fprintf (helpout, "    Output anchor segments, for reprocessing with --segments=<file>.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "rdotplot\n");
	fprintf (helpout, "    R output creates a file that can be plotted in the statistical package R.\n");
	fprintf (helpout, "    After creating the file like this:\n");
	fprintf (helpout, "        lastz ... --format=rdotplot > rdots.dat\n");
	fprintf (helpout, "    ask R to plot it using an R command like this:\n");
	fprintf (helpout, "        plot(read.table(\"rdots.dat\",header=T),type=\"l\")\n");
	fprintf (helpout, "    The separate option --rdotplot=<file> can be used to create a dot plot file\n");
	fprintf (helpout, "    at the same time as creating alignment output in another format.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "text\n");
	fprintf (helpout, "    Textual output is intended to be human readable.  Each alignment block is\n");
	fprintf (helpout, "    displayed with gap characters and a row of match/transition characters.\n");
	fprintf (helpout, "    Lines are wrapped at some reasonable width to allow printing to paper.\n");
	fprintf (helpout, "    The exact format of textual output may change in future releases of lastz.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "general\n");
	fprintf (helpout, "    General output creates a tab-delimited table with one line per alignment\n");
	fprintf (helpout, "    block.  The user can specify which fields are written (and in what order).\n");
	fprintf (helpout, "    This format is well-suited for use with spreadsheets and the R statistical\n");
	fprintf (helpout, "    package.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "    The format of the general output option is one of these:\n");
	fprintf (helpout, "        --format=general\n");
	fprintf (helpout, "        --format=general:fields\n");
	fprintf (helpout, "    where fields is a comma-separated list of field names.  If this list is\n");
	fprintf (helpout, "    absent all fields are printed.  The recognized field names are shown below.\n");
	fprintf (helpout, "    See the lastz readme file for more details.\n");
	fprintf (helpout, "\n");
	fprintf (helpout, "    The option --format=general- (with or without fields) inhibits the header\n");
	fprintf (helpout, "    lines.  This makes it suitable for catenating output from multiple runs.\n");

	lineWidth = 0;
	for (ix=0 ; genpafName[ix].name!=NULL ; ix++)
		{
		name    = genpafName[ix].name;
		nameLen = strlen(name);

		if (strcmp (name, "NA") == 0) continue;
		if (strcmp (name, "~")  == 0) continue;

		needComma = true;
		if (lineWidth == 0)
			{
			fprintf (helpout, "        ");
			lineWidth = 8;
			needComma = false;
			}
		else if (lineWidth + 2 + nameLen >= 79)
			{
			fprintf (helpout, ",\n        ");
			lineWidth = 8;
			needComma = false;
			}

		if (needComma) fprintf (helpout, ", ");
		fprintf (helpout, "%s", name);
		if (needComma) lineWidth += 2;
		lineWidth += nameLen;
		}
	fprintf (helpout, "\n");

	fprintf (helpout, "\n");
	fprintf (helpout, "The option --markend can be useful in cases (such as batch servers) in which\n");
	fprintf (helpout, "there may be a question as to whether or not lastz completed successfully.  The\n");
	fprintf (helpout, "line \"# lastz end-of-file\" is written to output as the last line.  Note that\n");
	fprintf (helpout, "in some formats this is *not* a legal line;  the user must remove it before any\n");
	fprintf (helpout, "downstream processsing.\n");

	exit (EXIT_FAILURE);
	}

static void shortcuts (void);
static void shortcuts (void)
	{
	scoreset* hoxScoring = NULL;

	fprintf (helpout, "%54s%s\n", "", "[defaults]");
	fprintf (helpout, "  B=0          same as --strand=plus%s\n",
	                 (defaultParams.whichStrand == 0)? "                   [B=0]" : "");
	fprintf (helpout, "  B=2          same as --strand=both\%s\n",
	                 (defaultParams.whichStrand > 0)? "                  [B=2]" : "");
	fprintf (helpout, "  B=-1         same as --strand=minus%s\n",
	                 (defaultParams.whichStrand < -1)? "                  [B=-1]" : "");

	fprintf (helpout, "  C=0          same as --nochain --gapped%s\n",
	                 ((!defaultParams.chain) && (defaultParams.gappedExtend))? "             [C=0]" : "");
	fprintf (helpout, "  C=1          same as --chain   --nogapped%s\n",
	                 ((defaultParams.chain) && (!defaultParams.gappedExtend))? "           [C=1]" : "");
	fprintf (helpout, "  C=2          same as --chain   --gapped%s\n",
	                 ((defaultParams.chain) && (defaultParams.gappedExtend))? "             [C=2]" : "");
	fprintf (helpout, "  C=3          same as --nochain --nogapped%s\n",
	                 ((!defaultParams.chain) && (!defaultParams.gappedExtend))? "           [C=3]" : "");

	fprintf (helpout, "  c=1          same as --census%23s[c=%d]\n",
	                 "", (defaultParams.reportCensus)? 1 : 0);

	fprintf (helpout, "  E=<penalty>  same as --gap=<..,penalty>%13s[E=" scoreFmtSimple "]\n",
	                 "", HOXD70_extend);

	fprintf (helpout, "  G=<score>    same as --chain=<score,..>%13s[G=" scoreFmtSimple "]\n",
	                 "", defaultParams.chainDiag);

	fprintf (helpout, "  H=<score>    same as --inner=<score>%16s[H=" scoreFmtSimple "]\n",
	                 "", defaultParams.innerThreshold);

	fprintf (helpout, "  K=<score>    same as --hspthresh=<score>%12s[K=%s]\n",
	                 "", score_thresh_to_string (&defaultParams.hspThreshold));
	fprintf (helpout, "  L=<score>    same as --gappedthresh=<score>%9s[L=K]\n", "");

	fprintf (helpout, "  M=<count>    same as --masking=<count>%14s[M=%d]\n",
	                 "", defaultParams.dynamicMasking);

	fprintf (helpout, "  m=<bytes>    same as --traceback=<bytes>%12s[m=%s]\n",
	                 "", unitize(defaultParams.tracebackMem));

	fprintf (helpout, "  O=<penalty>  same as --gap=<penalty,..>%13s[O=" scoreFmtSimple "]\n",
	                 "", HOXD70_open);

	fprintf (helpout, "  P=0          same as --noentropy%s\n",
	                 (!defaultParams.entropicHsp)? "                    [P=0]" : "");
	fprintf (helpout, "  P=1          same as --entropy%s\n",
	                 ((defaultParams.entropicHsp) && (!defaultParams.reportEntropy))? "                      [P=1]" : "");
	fprintf (helpout, "  P>1          same as --entropy=report%s\n",
	                 ((defaultParams.entropicHsp) && (defaultParams.reportEntropy))? "               [P>1]" : "");

	fprintf (helpout, "  Q=<file>     same as --scores=<file>%16s[Q=<HOXD70>]\n","");

	fprintf (helpout, "  R=<score>    same as --chain=<..,score>%13s[R=" scoreFmtSimple "]\n",
	                 "", defaultParams.chainAnti);

	fprintf (helpout, "  T=1          same as --seed=12of19 --transition%s\n",
	                 ((defaultParams.withTrans == 1) && (strcmp(defaultSeedString,seed_12of19) == 0))? "     [T=1]" : "");
	fprintf (helpout, "  T=2          same as --seed=12of19 --notransition%s\n",
	                 ((defaultParams.withTrans == 0) && (strcmp(defaultSeedString,seed_12of19) == 0))? "   [T=2]" : "");
	fprintf (helpout, "  T=3          same as --seed=14of22 --transition%s\n",
	                 ((defaultParams.withTrans == 1) && (strcmp(defaultSeedString,seed_14of22) == 0))? "     [T=3]" : "");
	fprintf (helpout, "  T=4          same as --seed=14of22 --notransition%s\n",
	                 ((defaultParams.withTrans == 0) && (strcmp(defaultSeedString,seed_14of22) == 0))? "   [T=4]" : "");

	fprintf (helpout, "  U=1          same as --match=1,1\n");

	fprintf (helpout, "  W=<length>   same as --seed=match<length>\n");

	fprintf (helpout, "  X=<score>    same as --xdrop=<score>%16s[X=10sub[A][A]]\n", "");
	fprintf (helpout, "  Y=<score>    same as --ydrop=<score>%16s[Y=O+300E]\n", "");

	fprintf (helpout, "  Z=<length>   same as --step=<length>%16s[Z=%u]\n",
	                 "", defaultParams.step);

	fprintf (helpout, "  v=0          same as --verbosity=0%s\n",
	                 (defaultParams.verbosity == 0)? "                  [v=0]" : "");
	fprintf (helpout, "  v=1          same as --verbosity=10%s\n",
	                 (defaultParams.verbosity == 10)? "                 [v=1]" : "");

	fprintf (helpout, "<HOXD70>\n");
	hoxScoring = new_dna_score_set (HOXD70, 0, 0, 0, 0);
	print_score_matrix_lf (helpout, hoxScoring, false, '\n');
	free_score_set ("<HOXD70>", hoxScoring);

	exit (EXIT_FAILURE);
	}


static void expander_options (char* header, char* prefix);
static void expander_options (char* header, char* prefix)
	{
	int ix, width, len;

	width = 0;
	for (ix=0 ; ix<numExpanders ; ix++)
		{
		if (strcmp_prefix (expanders[ix].argName, prefix) != 0) continue;
		len = strlen (expanders[ix].argName);
		if (len > width) width = len;
		}

	fprintf (helpout, "%s\n", header);
	if (width == 0)
		{
		fprintf (helpout, "  (none)\n");
		exit (EXIT_FAILURE);
		}

	for (ix=0 ; ix<numExpanders ; ix++)
		{
		if (strcmp_prefix (expanders[ix].argName, prefix) != 0) continue;
		fprintf (helpout, "  %-*s  (%s)\n",
		                  width, expanders[ix].argName,
		                  expanders[ix].expansion);
		}

	exit (EXIT_FAILURE);
	}


static void chastise (const char* format, ...);
static void chastise (const char* format, ...)
	{
	va_list	args;

	va_start (args, format);
	if (format != NULL)
		vfprintf (stderr, format, args);
	va_end (args);

	usage ();
	}


// 'local' variables shared by parse_options() and parse_options_loop()

char*	seq1Actions;
char*	seq2Actions;
char*	seedString;
char*	scoreFilename;
char*	infControlFilename;
int		haveXDrop;
int		haveYDrop;
int		haveStep;
int		haveGappedOption;
int		haveHspThreshold;
int		haveGappedThreshold;
int		haveInterpThreshold;
int		haveEntropicHsp;
int		haveBallScore;
int		haveWithTrans, haveWithTransForMatch;
int		haveMaxIdentity;
int		useUnitScores;
int		gappedExtendVerbosity;
int		unitMatch;
int		unitMismatch;
int		haveGapOpen;
int		haveGapExtend;
char*	gapOpenStr;
char*	gapExtendStr;
score	gapOpen;
score	gapExtend;
int		twinsYes, minGap, maxGap;
double	ballScoreFactor;
char*	firstSpecialSub;
score	specialSubScores[4][4];
int		formatIsSegments;
int		formatIsDotPlot;


// parse_options_loop-- main options parsing loop;  it is separate from the
//   overall options parsing routine because some options expand to a second
//   level of parsing

static void parse_options_loop   (int argc, char** argv,
                                  control* lzParams, control* izParams,
                                  int isToplevel);
static void parse_options_string (char* s,
                                  control* lzParams, control* izParams);

static void parse_options_loop
   (int			argc,
	char**		argv,
	control*	lzParams,
	control*	izParams,
	int			isTopLevel)
	{
	char*		argStr;
	int			argLen, argsLen;
	char*		scan, *scan2;
	int			items, scanned;
	char		extra;
	int			wordLen, ix;
	float		minPctId,      maxPctId;
	float		minCov,        maxCov;
	float		minContinuity, maxContinuity;
	int			r, c;
	int			argIsAMatch, charsUsed, tempInt;
	char*		waywardBracketArg = NULL;

	argsLen = 0; // (placate compiler)

	while (argc > 0)
		{
		// copy argument (if it turns out to be a file name, we'll erase it
		// later)

		if (isTopLevel)
			{
			argsLen = strlen(lzParams->args);
			strcpy (lzParams->args+argsLen, argv[0]);
			strcpy (lzParams->args+argsLen+strlen(argv[0])," ");
			}

		// locate arg value string, if there is one

		argStr = strchr(argv[0],'=');
		if (argStr != NULL) argStr++;

		// --svn (unadvertised)
		// $$$ This needs to be improved so that it shows the *latest* revision
		// $$$ .. number of any module in the build path

		if (strcmp (argv[0], "--svn") == 0)
			{
			printf ("SVN revision: %s\n", svnRevisionNumber);
			exit (EXIT_FAILURE);
			}

		// --self

		if (strcmp (argv[0], "--self") == 0)
			{
			lzParams->selfCompare = lzParams->inhibitTrivial = true;
			goto next_arg;
			}

		// --notrivial

		if (strcmp (argv[0], "--notrivial") == 0)
			{ lzParams->inhibitTrivial = true;  goto next_arg; }

		// --seed=<seed>, --seed=match<length> and variants

		if (strcmp (argv[0], "T=0") == 0)	// in blastz, T=0 was accompanied
			{								// .. by W=<length>, which did not
			lzParams->withTrans = 0;		// .. support transitions
			goto next_arg;
			}

		if (strcmp (argv[0], "T=1") == 0)
			{
			if (seedString != NULL) goto duplicated_option;
			seedString = copy_string (seed_12of19);
			lzParams->withTrans = 1;
			goto next_arg;
			}

		if (strcmp (argv[0], "T=2") == 0)
			{
			if (seedString != NULL) goto duplicated_option;
			seedString = copy_string (seed_12of19);
			lzParams->withTrans = 0;
			goto next_arg;
			}

		if (strcmp (argv[0], "T=3") == 0)
			{
			if (seedString != NULL) goto duplicated_option;
			seedString = copy_string (seed_14of22);
			lzParams->withTrans = 1;
			goto next_arg;
			}

		if (strcmp (argv[0], "T=4") == 0)
			{
			if (seedString != NULL) goto duplicated_option;
			seedString = copy_string (seed_14of22);
			lzParams->withTrans = 0;
			goto next_arg;
			}

		if (strcmp_prefix (argv[0], "W=") == 0)
			{
			if (seedString != NULL)
				chastise ("can't specify W=<length> with --seed\n");
			items = sscanf (argStr, "%d", &wordLen);
			if (items != 1) goto cant_understand;
			goto build_match_seed;
			}

		if (strcmp_prefix (argv[0], "--seed=") == 0)
			{
			if (seedString != NULL) goto duplicated_option;
			if (strcmp (argStr, "12of19") == 0)
				{
				seedString = copy_string (seed_12of19);
				goto next_arg;
				}

			if (strcmp (argStr, "14of22") == 0)
				{
				seedString = copy_string (seed_14of22);
				goto next_arg;
				}

			if ((strcmp_prefix (argStr, "match(") == 0)
			 && (argStr[strlen(argStr)-1] == ')'))
				{
				scan  = strchr(argStr,'(') + 1;
				items = sscanf (scan, "%d)", &wordLen);
				if (items != 1) goto cant_understand;
				goto build_match_seed;
				}

			if (strcmp_prefix (argStr, "match") == 0)
				{
				scan  = argStr + strlen("match");
				items = sscanf (scan, "%d%c", &wordLen, &extra);
				if (items != 1) goto cant_understand;
			build_match_seed:
				if ((wordLen < 2) || (wordLen > 15))
					chastise ("%d is not a valid word length\n", wordLen);

				seedString = malloc_or_die ("parse_options_loop (wordLen)", wordLen + 1);
				for (ix=0 ; ix<wordLen ; ix++)
					seedString[ix] = '1';
				seedString[wordLen] = 0;
				if (!haveWithTrans)
					{
					lzParams->withTrans = 0;
					haveWithTrans = true;
					haveWithTransForMatch = true;
					}
				goto next_arg;
				}

			if ((strcmp_prefix (argStr, "half(") == 0)
			 && (argStr[strlen(argStr)-1] == ')'))
				{
				scan  = strchr(argStr,'(') + 1;
				items = sscanf (scan, "%d)", &wordLen);
				if (items != 1) goto cant_understand;
				goto build_half_seed;
				}

			if (strcmp_prefix (argStr, "half") == 0)
				{
				scan  = argStr + strlen("half");
				items = sscanf (scan, "%d%c", &wordLen, &extra);
				if (items != 1) goto cant_understand;
			build_half_seed:
				if ((wordLen < 2) || (wordLen > 31))
					chastise ("%d is not a valid word length\n", wordLen);

				seedString = malloc_or_die ("parse_options_loop (wordLen)", wordLen + 1);
				for (ix=0 ; ix<wordLen ; ix++)
					seedString[ix] = 'T';
				seedString[wordLen] = 0;
				goto next_arg;
				}

			seedString = copy_string (argStr);
			goto next_arg;
			}

		// --[no]trans[ition]

		if ((strcmp (argv[0], "--notrans")      == 0)
		 || (strcmp (argv[0], "--notransition") == 0)
		 || (strcmp (argv[0], "--trans=0")      == 0)
		 || (strcmp (argv[0], "--transition=0") == 0))
			{ lzParams->withTrans = 0;  goto next_arg; }

		if ((strcmp (argv[0], "--trans")        == 0)
		 || (strcmp (argv[0], "--transition")   == 0)
		 || (strcmp (argv[0], "--trans=1")      == 0)
		 || (strcmp (argv[0], "--transition=1") == 0))
			{ lzParams->withTrans = 1;  haveWithTrans = true;  goto next_arg; }

		if ((strcmp (argv[0], "--trans=2")       == 0)
		 || (strcmp (argv[0], "--transition=2")  == 0)
		 || (strcmp (argv[0], "--transitions=2") == 0))
			{ lzParams->withTrans = 2;  haveWithTrans = true;  goto next_arg; }

		// --[no]filter=<[T,]M> and --filter=cares:<[T,]M>
		// (--filter=<[T:]M> supported for historical reasons)

		if (strcmp (argv[0], "--nofilter") == 0)
			{ lzParams->minMatches = -1;  goto next_arg; }

		if (strcmp_prefix (argv[0], "--filter=cares:") == 0)
			{
			scan = strchr(argStr,',');
			if (scan != NULL)
				{
				*(scan++) = 0;
				lzParams->maxTransversions = string_to_int (argStr + strlen("cares:"));
				lzParams->minMatches       = string_to_int (scan);
				}
			else
				{
				lzParams->maxTransversions = -1;
				lzParams->minMatches       = string_to_int (argStr + strlen("cares:"));
				}

			lzParams->filterCaresOnly = true;
			goto next_arg;
			}

		if (strcmp_prefix (argv[0], "--filter=") == 0)
			{
			scan = strchr(argStr,',');
			if (scan != NULL)
				{
				*(scan++) = 0;
				lzParams->maxTransversions = string_to_int (argStr);
				lzParams->minMatches       = string_to_int (scan);
				}
			else if ((scan = strchr(argStr,':')) != NULL)
				{
				*(scan++) = 0;
				lzParams->maxTransversions = string_to_int (argStr);
				lzParams->minMatches       = string_to_int (scan);
				}
			else
				{
				lzParams->maxTransversions = -1;
				lzParams->minMatches       = string_to_int (argStr);
				}

			lzParams->filterCaresOnly = false;
			goto next_arg;
			}

		// --word=<bits>

		if (strcmp_prefix (argv[0], "--word=") == 0)
			{
			lzParams->maxIndexBits = string_to_int (argStr);
			goto next_arg;
			}

		// --notwins and --twins=<min..max>
		// (--twins=<min:max> supported for historical reasons)

		if (strcmp (argv[0], "--notwins") == 0)
			{ twinsYes = false;  goto next_arg; }

		if (strcmp_prefix (argv[0], "--twins=") == 0)
			{
			twinsYes = true;
			scan = strstr(argStr,"..");
			if (scan != NULL)
				{
				*scan = 0;  scan += 2;
				minGap = string_to_int (argStr);
				maxGap = string_to_int (scan);
				}
			else if ((scan = strchr(argStr,':')) != NULL)
				{
				*(scan++) = 0;
				minGap = string_to_int (argStr);
				maxGap = string_to_int (scan);
				}
			else
				{
				minGap = 0;
				maxGap = string_to_int (argStr);
				}

			goto next_arg;
			}

		// --seedqueue=<entries>

#ifndef noSeedHitQueue
		if (strcmp_prefix (argv[0], "--seedqueue=") == 0)
			{
			lzParams->seedHitQueueSize = string_to_int (argStr);
			goto next_arg;
			}
#endif // not noSeedHitQueue

		// --recoverseeds (used to be --recoverhits)

		if ((strcmp (argv[0], "--norecoverseeds") == 0)
		 || (strcmp (argv[0], "--norecoverhits")  == 0))
			{ lzParams->basicHitType = hitSimple;  goto next_arg; }

		if ((strcmp (argv[0], "--recoverseeds") == 0)
		 || (strcmp (argv[0], "--recoverhits")  == 0))
			{ lzParams->basicHitType = hitRecover;  goto next_arg; }

		// --rawhits

		if (strcmp (argv[0], "--rawhits") == 0)
			{ lzParams->noHitFiltering = true;  goto next_arg; }

		// --step=<length> or Z=<length>

		if ((strcmp_prefix (argv[0], "--step=") == 0)
		 || (strcmp_prefix (argv[0], "Z=")      == 0))
			{
			tempInt = string_to_int (argStr);
			if (tempInt <= 0)
				suicidef ("--step must be positive");
			lzParams->step = tempInt;
			haveStep       = true;
			goto next_arg;
			}

		// --strand=both, etc.

		if ((strcmp (argv[0], "--both")        == 0)
		 || (strcmp (argv[0], "--bothstrands") == 0)
		 || (strcmp (argv[0], "--strand=both") == 0))
			{ lzParams->whichStrand = 2;  goto next_arg; }

		if ((strcmp (argv[0], "--plus")        == 0)
		 || (strcmp (argv[0], "--plusstrand")  == 0)
		 || (strcmp (argv[0], "--strand=plus") == 0))
			{ lzParams->whichStrand = 0;  goto next_arg; }

		if ((strcmp (argv[0], "--minus")        == 0)
		 || (strcmp (argv[0], "--minusstrand")  == 0)
		 || (strcmp (argv[0], "--strand=minus") == 0))
			{ lzParams->whichStrand = -1;  goto next_arg; }

		if (strcmp_prefix (argv[0], "B=") == 0)
			{
			lzParams->whichStrand = string_to_int (argStr);
			goto next_arg;
			}

		// --ambiguous=n and --ambiguous=iupac

		if ((strcmp (argv[0], "--ambiguousn")  == 0)
		 || (strcmp (argv[0], "--ambig=n")     == 0)
		 || (strcmp (argv[0], "--ambiguous=n") == 0))
			{ lzParams->nIsAmbiguous = true;  goto next_arg; }

		if ((strcmp (argv[0], "--ambiguous=iupac") == 0)
		 || (strcmp (argv[0], "--ambig=iupac")     == 0))
			{
			lzParams->allowAmbiDNA = lzParams->nIsAmbiguous = true;
			goto next_arg;
			}

		// --[no]gfextend or (unadvertised) --[no]gfx

		if ((strcmp (argv[0], "--gfextend") == 0)
		 || (strcmp (argv[0], "--gfx"     ) == 0))
			{ lzParams->gfExtend = gfexXDrop;  goto next_arg; }

		if ((strcmp (argv[0], "--nogfextend") == 0)
		 || (strcmp (argv[0], "--nogfx"     ) == 0))
			{ lzParams->gfExtend = gfexNoExtend;  goto next_arg; }

		// --justhits or --hitsonly (unadvertised)

		if ((strcmp (argv[0], "--justhits") == 0)
		 || (strcmp (argv[0], "--hitsonly") == 0))
			{
			lzParams->gfExtend     = gfexNoExtend;
			lzParams->gappedExtend = false;
			goto next_arg;
			}

		// --[no]chain, --chain=<diag,anti>, G=, or R=

		if (strcmp (argv[0], "--chain") == 0)
			{ lzParams->chain = true;  goto next_arg; }

		if (strcmp (argv[0], "--nochain") == 0)
			{ lzParams->chain = false;  goto next_arg; }

		if (strcmp_prefix (argv[0], "--chain=") == 0)
			{
			lzParams->chain = true;
			scan = strchr(argStr,',');
			if (scan == NULL)
				chastise ("%s is not a valid pair of chain penalties\n", argStr);
			*scan = 0;
			lzParams->chainDiag = string_to_score (argStr);
			*(scan++) = ',';
			lzParams->chainAnti = string_to_score (scan);
			goto next_arg;
			}

		if (strcmp_prefix (argv[0], "G=") == 0)
			{
			lzParams->chainDiag = string_to_score (argStr);
			goto next_arg;
			}

		if (strcmp_prefix (argv[0], "R=") == 0)
			{
			lzParams->chainAnti = string_to_score (argStr);
			goto next_arg;
			}

		// --[no]gapped, C=<code>, or (unadvertised) --[no]gx

		if ((strcmp (argv[0], "--gapped") == 0)
		 || (strcmp (argv[0], "--gx"    ) == 0))
			{
			lzParams->gappedExtend = true;
			haveGappedOption = true;
			goto next_arg;
			}

		if ((strcmp (argv[0], "--nogapped") == 0)
		 || (strcmp (argv[0], "--ungapped") == 0)
		 || (strcmp (argv[0], "--nogx"    ) == 0))
			{ lzParams->gappedExtend = false;  goto next_arg; }

		if (strcmp (argv[0], "C=0") == 0)
			{
			lzParams->chain = false;
			lzParams->gappedExtend = true;
			haveGappedOption = true;
			goto next_arg;
			}

		if (strcmp (argv[0], "C=1") == 0)
			{
			lzParams->chain = true;
			lzParams->gappedExtend = false;
			goto next_arg;
			}

		if (strcmp (argv[0], "C=2") == 0)
			{
			lzParams->chain = true;
			lzParams->gappedExtend = true;
			haveGappedOption = true;
			goto next_arg;
			}

		if (strcmp (argv[0], "C=3") == 0)
			{
			lzParams->chain = false;
			lzParams->gappedExtend = false;
			goto next_arg;
			}

		// --stopafterone

		if ((strcmp (argv[0], "--stopafterone") == 0)
		 || (strcmp (argv[0], "--anyornone")    == 0))
			{
			lzParams->hspToGapped = true;
			lzParams->searchLimit = 1;
			goto next_arg;
			}

		// --allgappedbounds

		if (strcmp (argv[0], "--allgappedbounds") == 0)
			{ lzParams->gappedAllBounds = true;  goto next_arg; }

		// --score[s]=<file> or Q=<file>

		if ((strcmp_prefix (argv[0], "--scores=") == 0)
		 || (strcmp_prefix (argv[0], "--score=")  == 0)
		 || (strcmp_prefix (argv[0], "Q=")        == 0))
			{
			if (scoreFilename != NULL) goto duplicated_option;
			scoreFilename = copy_string (argStr);
			goto next_arg;
			}

		// --match=<match>[,<mismatch>]
		// or (unadvertised) --unitscore[s] or U=1
		// or (unadvertised) U=<match>[,<mismatch>]

		if ((strcmp (argv[0], "--unitscore")  == 0)
		 || (strcmp (argv[0], "--unitscores") == 0)
		 || (strcmp (argv[0], "U=1")          == 0))
			{
			useUnitScores = true;
			unitMatch     = 1;
			unitMismatch  = -1;
			goto next_arg;
			}

		if ((strcmp_prefix (argv[0], "--match=") == 0)
		 || (strcmp_prefix (argv[0], "U=")       == 0))
			{
			useUnitScores = true;
			scan = strchr(argStr,',');
			if (scan != NULL)
				*(scan++) = 0;
			unitMatch = string_to_score (argStr);
			if (unitMatch <= 0)
				chastise ("%s is not a valid match score\n", argStr);
			if (scan == NULL)
				unitMismatch = -unitMatch;
			else
				{
				unitMismatch = -string_to_score (scan);
				if (unitMismatch >= 0)
					chastise ("%s is not a valid mismatch penalty\n", scan);
				}
			goto next_arg;
			}

		// (unadvertised) <nuc1><nuc2>:<score>

		if ((nuc_to_bits[((u8*)argv[0])[0]] >= 0)
		 && (nuc_to_bits[((u8*)argv[0])[1]] >= 0)
		 && (argv[0][2] == ':'))
			{
			r = nuc_to_bits[(u8)argv[0][0]];
			c = nuc_to_bits[(u8)argv[0][1]];
			specialSubScores[r][c] = string_to_score (&argv[0][3]);
			if (firstSpecialSub == NULL)
				firstSpecialSub = copy_string (argv[0]);
			goto next_arg;
			}

		// --infer[only][=<control>]

		if (strcmp (argv[0], "--infer") == 0)
			{
			if (infControlFilename != NULL) goto duplicated_option;
			lzParams->inferScores = true;
			lzParams->inferOnly   = false;
			goto next_arg;
			}

		if (strcmp_prefix (argv[0], "--infer=") == 0)
			{
			if (infControlFilename != NULL) goto duplicated_option;
			lzParams->inferScores = true;
			lzParams->inferOnly   = false;
			infControlFilename = copy_string (argStr);
			goto next_arg;
			}

		if (strcmp (argv[0], "--inferonly") == 0)
			{
			if (infControlFilename != NULL) goto duplicated_option;
			lzParams->inferScores = true;
			lzParams->inferOnly   = true;
			goto next_arg;
			}

		if (strcmp_prefix (argv[0], "--inferonly=") == 0)
			{
			if (infControlFilename != NULL) goto duplicated_option;
			lzParams->inferScores = true;
			lzParams->inferOnly   = true;
			infControlFilename = copy_string (argStr);
			goto next_arg;
			}

		// --infscores[=<file>]

		if (strcmp (argv[0], "--infscores") == 0)
			{
			lzParams->inferScores = true;
			goto next_arg;
			}

		if (strcmp_prefix (argv[0], "--infscores=") == 0)
			{
			if (izParams->ic.inferFilename != NULL) goto duplicated_option;
			lzParams->inferScores      = true;
			izParams->ic.inferFilename = copy_string (argStr);
			goto next_arg;
			}

		// --gap=<[open,]extend> or O=<open> or E=<extend>

		if (strcmp_prefix (argv[0], "--gap=") == 0)
			{
			scan = strchr(argStr,',');
			if (scan == NULL)
				{
				gapExtend = string_to_score (argStr);
				if (gapExtend < 0)
					chastise ("%s is not a valid gap extension penalty\n", argStr);
				haveGapExtend = true;
				gapExtendStr  = argStr;
				}
			else
				{
				*(scan++) = 0;
				gapOpen      = string_to_score (argStr);
				gapExtend    = string_to_score (scan);
				haveGapOpen  = haveGapExtend = true;
				gapOpenStr   = argStr;
				gapExtendStr = scan;
				}
			goto next_arg;
			}

		if (strcmp_prefix (argv[0], "O=") == 0)
			{
			gapOpen = string_to_score (argStr);
			haveGapOpen = true;
			gapOpenStr  = argStr;
			goto next_arg;
			}

		if (strcmp_prefix (argv[0], "E=") == 0)
			{
			gapExtend = string_to_score (argStr);
			haveGapExtend = true;
			gapExtendStr  = argStr;
			goto next_arg;
			}

		// --xdrop=<score> or X=<score>

		if ((strcmp_prefix (argv[0], "--xdrop=") == 0)
		 || (strcmp_prefix (argv[0], "X=")       == 0))
			{
			if ((haveHspThreshold)
			 && (lzParams->gfExtend == gfexExact))
				chastise ("can't use %s with --exact\n", argv[0]);
			if ((haveHspThreshold)
			 && (lzParams->gfExtend >= gfexMismatch_min)
			 && (lzParams->gfExtend <= gfexMismatch_max))
				chastise ("can't use %s with --%dmismatch\n", argv[0], lzParams->gfExtend);
			lzParams->gfExtend = gfexXDrop;
			lzParams->xDrop    = string_to_score (argStr);
			haveXDrop = true;
			goto next_arg;
			}

		// --ydrop=<score> or Y=<score>

		if ((strcmp_prefix (argv[0], "--ydrop=") == 0)
		 || (strcmp_prefix (argv[0], "Y=")       == 0))
			{
			lzParams->yDrop = string_to_score (argStr);
			haveYDrop = true;
			goto next_arg;
			}

		// --noxtrim

		if ((strcmp (argv[0], "--noxtrim")     == 0)
		 || (strcmp (argv[0], "--noxdroptrim") == 0))
			{ chastise ("sorry, --noxtrim not implemented yet\n", argv[0]);
			  lzParams->xDropUntrimmed = true;  goto next_arg; }

		// --noytrim

		if ((strcmp (argv[0], "--noytrim")     == 0)
		 || (strcmp (argv[0], "--noydroptrim") == 0))
			{ lzParams->yDropUntrimmed = true;  goto next_arg; }

		// --hspthresh=<score_threshold> or K=<score_threshold>

		if ((strcmp_prefix (argv[0], "--hspthresh=")    == 0)
		 || (strcmp_prefix (argv[0], "--hspthreshold=") == 0)
		 || (strcmp_prefix (argv[0], "--mspthresh=")    == 0)
		 || (strcmp_prefix (argv[0], "--mspthreshold=") == 0)
		 || (strcmp_prefix (argv[0], "K="          )    == 0))
			{
			if ((haveHspThreshold)
			 && (lzParams->gfExtend == gfexExact))
				chastise ("can't use %s with --exact\n", argv[0]);
			if ((haveHspThreshold)
			 && (lzParams->gfExtend >= gfexMismatch_min)
			 && (lzParams->gfExtend <= gfexMismatch_max))
				chastise ("can't use %s with --%dmismatch\n", argv[0], lzParams->gfExtend);
			lzParams->hspThreshold = string_to_score_thresh (argStr);
			haveHspThreshold = true;
			goto next_arg;
			}

		// --exact=<length>

		if (strcmp_prefix (argv[0], "--exact=") == 0)
			{
			scan = argStr;
		parse_gfex_exact:
			if ((haveHspThreshold)
			 && (lzParams->gfExtend == gfexXDrop))
				chastise ("can't use %s with --hspthreshold\n", argv[0]);
			if ((haveXDrop)
			 && (lzParams->gfExtend == gfexXDrop))
				chastise ("can't use %s with --xdrop\n", argv[0]);
			if ((haveHspThreshold)
			 && (lzParams->gfExtend >= gfexMismatch_min)
			 && (lzParams->gfExtend <= gfexMismatch_max))
				chastise ("can't use %s with --%dmismatch\n", argv[0], lzParams->gfExtend);
			lzParams->gfExtend = gfexExact;
			lzParams->hspThreshold.t = 'S';
			lzParams->hspThreshold.s = string_to_score (scan);
			if (lzParams->hspThreshold.s <= 0)
				chastise ("%s is not a valid exact match threshold\n", scan);
			haveHspThreshold = true;
			goto next_arg;
			}

		// --<N>mismatch=<length> or  --mismatch=<M>,<length>

		argIsAMatch = false;
		if (sscanf (argv[0], "--%d%n", &tempInt, &charsUsed) == 1)
			argIsAMatch = (strcmp_prefix (argv[0]+charsUsed, "mismatch=") == 0);

		if (argIsAMatch)
			{ scan = argStr;  goto parse_gfex_mismatch; }

		if (strcmp_prefix (argv[0], "--mismatch=") == 0)
			{
			scan = strchr(argStr,',');
			if (scan == NULL)
				chastise ("--mismatch requires two values (count and length)\n");
			*(scan++) = 0;
			tempInt = string_to_score (argStr);
		parse_gfex_mismatch:
			if (tempInt == 0) goto parse_gfex_exact;
			if ((tempInt < gfexMismatch_min) || (tempInt > gfexMismatch_max))
				chastise ("%d is out of range for N-mismatch (valid range is %d..%d)\n",
				          tempInt, gfexMismatch_min, gfexMismatch_max);
			if ((haveHspThreshold)
			 && (lzParams->gfExtend == gfexXDrop))
				chastise ("can't use %s with --hspthreshold\n", argv[0]);
			if ((haveXDrop)
			 && (lzParams->gfExtend == gfexXDrop))
				chastise ("can't use %s with --xdrop\n", argv[0]);
			if ((haveHspThreshold)
			 && (lzParams->gfExtend == gfexExact))
				chastise ("can't use %s with --exact\n", argv[0]);
			lzParams->gfExtend = tempInt;
			lzParams->hspThreshold.t = 'S';
			lzParams->hspThreshold.s = string_to_score (scan);
			if (lzParams->hspThreshold.s < tempInt)
				chastise ("%s is not a valid exact %dmismatch threshold\n", scan, tempInt);
			haveHspThreshold = true;
			goto next_arg;
			}

		// --inner=<score> or H=<score>

		if ((strcmp_prefix (argv[0], "--inner=") == 0)
		 || (strcmp_prefix (argv[0], "H="      ) == 0))
			{
			lzParams->innerThreshold = string_to_score (argStr);
			haveInterpThreshold = true;
			goto next_arg;
			}

		// --gappedthresh=<score_threshold> or L=<score_threshold>

		if ((strcmp_prefix (argv[0], "--gappedthresh=")    == 0)
		 || (strcmp_prefix (argv[0], "--gappedthreshold=") == 0)
		 || (strcmp_prefix (argv[0], "L="             )    == 0))
			{
			lzParams->gappedThreshold = string_to_score_thresh (argStr);
			haveGappedThreshold = true;
			goto next_arg;
			}

		// --ball=<score>

		if (strcmp_prefix (argv[0], "--ball=") == 0)
			{
			argLen = strlen(argStr);
			if ((argLen > 0) && (argStr[argLen-1] == '%'))
				{
				lzParams->ballScore = 0;		// (just signals that --ball used)
				ballScoreFactor = pct_string_to_double (argStr);
				}
			else
				{
				lzParams->ballScore = string_to_score (argStr);
				haveBallScore = true;
				}
			goto next_arg;
			}

		// --[no]entropy or P=<flag>

		if ((strcmp (argv[0], "--entropy") == 0)
		 || (strcmp (argv[0], "P=1"      ) == 0))
			{ lzParams->entropicHsp = haveEntropicHsp = true;  goto next_arg; }

		if ((strcmp (argv[0], "--noentropy") == 0)
		 || (strcmp (argv[0], "P=0"        ) == 0))
			{ lzParams->entropicHsp = haveEntropicHsp = false;  goto next_arg; }

		if (strcmp_prefix (argv[0], "P=") == 0)
			{
			if (string_to_int (argStr) <= 0)
				chastise ("illegal value for P");
			goto report_entropy;
			}

		if (strcmp (argv[0], "--entropy=report") == 0)
			{
		report_entropy:
			lzParams->entropicHsp = lzParams->reportEntropy = haveEntropicHsp = true;
			goto next_arg;
			}

		// --traceback=<bytes> or m=<bytes>

		if ((strcmp_prefix (argv[0], "--traceback=") == 0)
		 || (strcmp_prefix (argv[0], "m="          ) == 0))
			{
			lzParams->tracebackMem = string_to_int2 (argStr, false);
			goto next_arg;
			}

		// --maxwordcount=<limit>[%][,<chasm>]

		if (strcmp_prefix (argv[0], "--maxwordcount=") == 0)
			{
			argLen = strlen(argStr);
			scan   = strchr(argStr,',');
			if (scan != NULL)
				{
				argLen = scan - argStr;
				*(scan++) = 0;
				tempInt = string_to_int (scan);
				if (tempInt < 1)
					suicidef ("--maxwordcount's max interval must be at least 1");
				lzParams->maxWordCountChasm = tempInt;
				}
			if ((argLen > 0) && (argStr[argLen-1] == '%'))
				{
				lzParams->wordCountKeep  = pct_string_to_double (argStr);
				lzParams->wordCountLimit = 0;
				if (lzParams->wordCountKeep < 0)
					suicidef ("--maxwordcount cannot be zero");
				else if (lzParams->wordCountKeep < 0)
					suicidef ("--maxwordcount cannot be negative");
				else if (lzParams->wordCountKeep == 1)
					suicidef ("--maxwordcount cannot be 100%");
				else if (lzParams->wordCountKeep >= 1)
					suicidef ("--maxwordcount cannot be more than 100%");
				}
			else
				{
				tempInt = string_to_int (argStr);
				if (tempInt < 1)
					suicidef ("--maxwordcount must be at least 1");
				lzParams->wordCountLimit = tempInt;
				lzParams->wordCountKeep  = 0.0;
				}
			goto next_arg;
			}

		// --masking=<count> or M=<count>

		if ((strcmp_prefix (argv[0], "--masking=") == 0)
		 || (strcmp_prefix (argv[0], "M="        ) == 0))
			{
			tempInt = string_to_int (argStr);
			if (tempInt < 0)
				suicidef ("--masking cannot be negative");
			else if (tempInt >= 255)  // (255 itself is not allowed)
				suicidef ("--masking must be less than 255");
			lzParams->dynamicMasking = tempInt;
			goto next_arg;
			}

		// --[no]census[=<file>] or c=<flag>

		if (strcmp_prefix (argv[0], "c=") == 0)
			{
			if (string_to_int (argStr) == 0) goto census_off;
			                            else goto census_on;
			}

		if (strcmp (argv[0], "--census") == 0)
			{
		census_on:
			lzParams->reportCensus = true;
			lzParams->censusKind   = 'B';
			goto next_arg;
			}

		if (strcmp (argv[0], "--nocensus") == 0)
			{
		census_off:
			lzParams->reportCensus = false;
			goto next_arg;
			}

		if (strcmp_prefix (argv[0], "--census=") == 0)
			{
			if (lzParams->censusFilename != NULL) goto duplicated_option;
			lzParams->censusKind = 'B';
			goto census_to_file;
			}

		if (strcmp_prefix (argv[0], "--census16=") == 0)
			{
			if (lzParams->censusFilename != NULL) goto duplicated_option;
			lzParams->censusKind = 'W';
			goto census_to_file;
			}

		if (strcmp_prefix (argv[0], "--census32=") == 0)
			{
			if (lzParams->censusFilename != NULL) goto duplicated_option;
			lzParams->censusKind     = 'L';
		census_to_file:
			lzParams->censusFilename = copy_string (argStr);
			lzParams->reportCensus   = true;
			goto next_arg;
			}

		// --identity=min[..max]

		if (strcmp_prefix (argv[0], "--identity=") == 0)
			{
			if (strcmp (argStr,"..") == 0) goto cant_understand;
			scan  = argStr;
			scan2 = strstr (scan, "..");
			if      (scan2 == NULL) { ;                           } // min
			else if (scan2 == scan) { scan2 += 2;    scan = NULL; } // ..max
			else if (scan2[2] == 0) { scan2 = NULL;               } // min..
			else                    { scan2 += 2;                 } // min..max

			minPctId = 0.0;
			if (scan != NULL)
				{
				scanned = -1;
				sscanf (scan, "%f%n", &minPctId, &scanned);
				if (scanned == -1) goto cant_understand;
				scan += scanned;
				if (*scan == '%') scan++;
				if (scan2 == NULL)
					{
					if ((*scan != 0) && (strcmp (scan, ".") != 0) && (strcmp (scan, "..") != 0))
						goto cant_understand;
					}
				else
					{
					if (*scan != '.') goto cant_understand;  scan++;
					if (scan != scan2) { if (*scan != '.') goto cant_understand;  scan++; }
					if (scan != scan2) goto cant_understand;
					}
				}

			maxPctId = 100.0;
			if (scan2 != NULL)
				{
				scanned = -1;
				sscanf (scan2, "%f%n", &maxPctId, &scanned);
				if (scanned == -1) goto cant_understand;
				scan2 += scanned;
				if (*scan2 == '%') scan2++;
				if (*scan2 != 0) goto cant_understand;
				}

			if ((minPctId < 0) || (maxPctId > 100) || (minPctId > maxPctId))
				goto cant_understand;

			lzParams->minIdentity = minPctId / 100.0;
			lzParams->maxIdentity = maxPctId / 100.0;
			haveMaxIdentity       = true;
			goto next_arg;
			}

		// --coverage=min[..max]

		if (strcmp_prefix (argv[0], "--coverage=") == 0)
			{
			if (strcmp (argStr,"..") == 0) goto cant_understand;
			scan  = argStr;
			scan2 = strstr (scan, "..");
			if      (scan2 == NULL) { ;                           } // min
			else if (scan2 == scan) { scan2 += 2;    scan = NULL; } // ..max
			else if (scan2[2] == 0) { scan2 = NULL;               } // min..
			else                    { scan2 += 2;                 } // min..max

			minCov = 0.0;
			if (scan != NULL)
				{
				scanned = -1;
				sscanf (scan, "%f%n", &minCov, &scanned);
				if (scanned == -1) goto cant_understand;
				scan += scanned;
				if (*scan == '%') scan++;
				if (scan2 == NULL)
					{
					if ((*scan != 0) && (strcmp (scan, ".") != 0) && (strcmp (scan, "..") != 0))
						goto cant_understand;
					}
				else
					{
					if (*scan != '.') goto cant_understand;  scan++;
					if (scan != scan2) { if (*scan != '.') goto cant_understand;  scan++; }
					if (scan != scan2) goto cant_understand;
					}
				}

			maxCov = 100.0;
			if (scan2 != NULL)
				{
				scanned = -1;
				sscanf (scan2, "%f%n", &maxCov, &scanned);
				if (scanned == -1) goto cant_understand;
				scan2 += scanned;
				if (*scan2 == '%') scan2++;
				if (*scan2 != 0) goto cant_understand;
				}

			if ((minCov < 0) || (maxCov > 100) || (minCov > maxCov))
				goto cant_understand;

			lzParams->minCoverage = minCov / 100.0;
			lzParams->maxCoverage = maxCov / 100.0;
			goto next_arg;
			}

		// --continuity=min[..max]

		if (strcmp_prefix (argv[0], "--continuity=") == 0)
			{
			if (strcmp (argStr,"..") == 0) goto cant_understand;
			scan  = argStr;
			scan2 = strstr (scan, "..");
			if      (scan2 == NULL) { ;                           } // min
			else if (scan2 == scan) { scan2 += 2;    scan = NULL; } // ..max
			else if (scan2[2] == 0) { scan2 = NULL;               } // min..
			else                    { scan2 += 2;                 } // min..max

			minContinuity = 0.0;
			if (scan != NULL)
				{
				scanned = -1;
				sscanf (scan, "%f%n", &minContinuity, &scanned);
				if (scanned == -1) goto cant_understand;
				scan += scanned;
				if (*scan == '%') scan++;
				if (scan2 == NULL)
					{
					if ((*scan != 0) && (strcmp (scan, ".") != 0) && (strcmp (scan, "..") != 0))
						goto cant_understand;
					}
				else
					{
					if (*scan != '.') goto cant_understand;  scan++;
					if (scan != scan2) { if (*scan != '.') goto cant_understand;  scan++; }
					if (scan != scan2) goto cant_understand;
					}
				}

			maxContinuity = 100.0;
			if (scan2 != NULL)
				{
				scanned = -1;
				sscanf (scan2, "%f%n", &maxContinuity, &scanned);
				if (scanned == -1) goto cant_understand;
				scan2 += scanned;
				if (*scan2 == '%') scan2++;
				if (*scan2 != 0) goto cant_understand;
				}

			if ((minContinuity < 0) || (maxContinuity > 100) || (minContinuity > maxContinuity))
				goto cant_understand;

			lzParams->minContinuity = minContinuity / 100.0;
			lzParams->maxContinuity = maxContinuity / 100.0;
			goto next_arg;
			}

		// --matchcount=<min>

		if (strcmp_prefix (argv[0], "--matchcount=") == 0)
			{
			tempInt = string_to_int (argStr);
			if (tempInt <= 0)
				suicidef ("--matchcount must be positive");
			lzParams->minMatchCount = tempInt;
			goto next_arg;
			}

		// --density=max

#ifdef densityFiltering
		if (strcmp_prefix (argv[0], "--density=") == 0)
			{
			lzParams->maxDensity = string_to_double (argStr);
			goto next_arg;
			}
#endif // densityFiltering

		// --[no]mirror

		if (strcmp (argv[0], "--mirror") == 0)
			{ lzParams->mirror = true;  goto next_arg; }

		if (strcmp (argv[0], "--nomirror") == 0)
			{ lzParams->mirror = false;  goto next_arg; }

		// --out[put]=<file>

		if ((strcmp_prefix (argv[0], "--out=")    == 0)
		 || (strcmp_prefix (argv[0], "--output=") == 0))
			{ 
			if (lzParams->outputFilename != NULL)
				free_if_valid ("output file name", lzParams->outputFilename);
			lzParams->outputFilename = copy_string (argStr);
			goto next_arg;
			}

		// --format=<format> (many variations)

		if ((strcmp (argv[0], "--format=gfa") == 0)
		 || (strcmp (argv[0], "--format=GFA") == 0)
		 || (strcmp (argv[0], "--gfa")        == 0)
		 || (strcmp (argv[0], "--GFA")        == 0))
			{ lzParams->outputFormat = fmtGfa;  goto next_arg; }

		if ((strcmp (argv[0], "--format=gfanoscore") == 0) // (unadvertised)
		 || (strcmp (argv[0], "--format=GFANOSCORE") == 0)
		 || (strcmp (argv[0], "--gfanoscore")        == 0)
		 || (strcmp (argv[0], "--GFANOSCORE")        == 0))
			{ lzParams->outputFormat = fmtGfaNoScore;  goto next_arg; }

		if ((strcmp (argv[0], "--format=lav") == 0)
		 || (strcmp (argv[0], "--format=LAV") == 0)
		 || (strcmp (argv[0], "--lav")        == 0)
		 || (strcmp (argv[0], "--LAV")        == 0))
			{ lzParams->outputFormat = fmtLav;  goto next_arg; }

		if ((strcmp (argv[0], "--format=lav+") == 0)
		 || (strcmp (argv[0], "--format=LAV+") == 0)
		 || (strcmp (argv[0], "--lav+")        == 0)
		 || (strcmp (argv[0], "--LAV+")        == 0))
			{ lzParams->outputFormat = fmtLavComment;  goto next_arg; }

		if ((strcmp (argv[0], "--format=lav+text") == 0)
		 || (strcmp (argv[0], "--format=LAV+text") == 0)
		 || (strcmp (argv[0], "--lav+text")        == 0)
		 || (strcmp (argv[0], "--LAV+text")        == 0)
		 || (strcmp (argv[0], "--format=text+lav") == 0)
		 || (strcmp (argv[0], "--format=text+LAV") == 0)
		 || (strcmp (argv[0], "--text+lav")        == 0)
		 || (strcmp (argv[0], "--text+LAV")        == 0))
			{ lzParams->outputFormat = fmtLavText;  goto next_arg; }

		if ((strcmp (argv[0], "--format=lavscore") == 0) // (unadvertised)
		 || (strcmp (argv[0], "--format=LAVSCORE") == 0)
		 || (strcmp (argv[0], "--lavscore")        == 0)
		 || (strcmp (argv[0], "--LAVSCORE")        == 0))
			{ lzParams->outputFormat = fmtLavScore;  goto next_arg; }

		if ((strcmp (argv[0], "--format=axt") == 0)
		 || (strcmp (argv[0], "--format=AXT") == 0)
		 || (strcmp (argv[0], "--axt")        == 0)
		 || (strcmp (argv[0], "--AXT")        == 0))
			{ lzParams->outputFormat = fmtAxt;  goto next_arg; }

		if ((strcmp (argv[0], "--format=axt+") == 0)
		 || (strcmp (argv[0], "--format=AXT+") == 0)
		 || (strcmp (argv[0], "--axt+")        == 0)
		 || (strcmp (argv[0], "--AXT+")        == 0))
			{ lzParams->outputFormat = fmtAxtComment;  goto next_arg; }

		if ((strcmp (argv[0], "--format=maf") == 0)
		 || (strcmp (argv[0], "--format=MAF") == 0)
		 || (strcmp (argv[0], "--maf")        == 0)
		 || (strcmp (argv[0], "--MAF")        == 0))
			{
			lzParams->outputFormat = fmtMaf;
			maf_distinguishNames   = false;
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=~maf") == 0)
		 || (strcmp (argv[0], "--format=~MAF") == 0))
			{
			lzParams->outputFormat = fmtMaf;
			maf_distinguishNames   = true;
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=maf+") == 0)
		 || (strcmp (argv[0], "--format=MAF+") == 0)
		 || (strcmp (argv[0], "--maf+")        == 0)
		 || (strcmp (argv[0], "--MAF+")        == 0))
			{
			lzParams->outputFormat = fmtMafComment;
			maf_distinguishNames   = false;
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=~maf+") == 0)
		 || (strcmp (argv[0], "--format=~MAF+") == 0))
			{
			lzParams->outputFormat = fmtMafComment;
			maf_distinguishNames   = true;
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=maf-") == 0)
		 || (strcmp (argv[0], "--format=MAF-") == 0)
		 || (strcmp (argv[0], "--maf-")        == 0)
		 || (strcmp (argv[0], "--MAF-")        == 0))
			{
			lzParams->outputFormat = fmtMafNoComment;
			maf_distinguishNames   = false;
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=mafsegments") == 0)
		 || (strcmp (argv[0], "--format=MAFSEGMENTS") == 0)
		 || (strcmp (argv[0], "--mafsegments")        == 0)
		 || (strcmp (argv[0], "--MAFSEGMENTS")        == 0))
			{
			lzParams->outputFormat   = fmtMaf;
			lzParams->deGapifyOutput = true;
			maf_distinguishNames     = false;
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=mafsegments+") == 0)
		 || (strcmp (argv[0], "--format=MAFSEGMENTS+") == 0)
		 || (strcmp (argv[0], "--mafsegments+")        == 0)
		 || (strcmp (argv[0], "--MAFSEGMENTS+")        == 0))
			{
			lzParams->outputFormat   = fmtMafComment;
			lzParams->deGapifyOutput = true;
			maf_distinguishNames     = false;
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=mafsegments-") == 0)
		 || (strcmp (argv[0], "--format=MAFSEGMENTS-") == 0)
		 || (strcmp (argv[0], "--mafsegments-")        == 0)
		 || (strcmp (argv[0], "--MAFSEGMENTS-")        == 0))
			{
			lzParams->outputFormat   = fmtMafNoComment;
			lzParams->deGapifyOutput = true;
			maf_distinguishNames     = false;
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=softsam") == 0)
		 || (strcmp (argv[0], "--format=SOFTSAM") == 0)
		 || (strcmp (argv[0], "--softsam")        == 0)
		 || (strcmp (argv[0], "--SOFTSAM")        == 0))
			{
			lzParams->outputFormat = fmtSoftSam;
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=softsam-") == 0)
		 || (strcmp (argv[0], "--format=SOFTSAM-") == 0)
		 || (strcmp (argv[0], "--softsam-")        == 0)
		 || (strcmp (argv[0], "--SOFTSAM-")        == 0))
			{
			lzParams->outputFormat = fmtSoftSamNoHeader;
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=sam") == 0)
		 || (strcmp (argv[0], "--format=SAM") == 0)
		 || (strcmp (argv[0], "--sam")        == 0)
		 || (strcmp (argv[0], "--SAM")        == 0))
			{
			lzParams->outputFormat = fmtHardSam;
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=sam-") == 0)
		 || (strcmp (argv[0], "--format=SAM-") == 0)
		 || (strcmp (argv[0], "--sam-")        == 0)
		 || (strcmp (argv[0], "--SAM-")        == 0))
			{
			lzParams->outputFormat = fmtHardSamNoHeader;
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=cigar") == 0)
		 || (strcmp (argv[0], "--format=CIGAR") == 0)
		 || (strcmp (argv[0], "--cigar")        == 0)
		 || (strcmp (argv[0], "--CIGAR")        == 0))
			{
			lzParams->outputFormat = fmtCigar;
			goto next_arg;
			}

		if (strcmp_prefix (argv[0], "--writesegments=") == 0)
			{
			if (lzParams->outputFilename != NULL)
				free_if_valid ("output file name", lzParams->outputFilename);
			lzParams->outputFilename = copy_string (argStr);
			goto format_segments;
			}

		if (strcmp (argv[0], "--format=segments") == 0) // (now unadvertised!)
			{
		format_segments:
			free_if_valid ("lzParams->outputInfo", lzParams->outputInfo);
			lzParams->outputFormat = fmtGenpaf;
			lzParams->outputInfo   = copy_string (genpafSegmentKeys);
			formatIsSegments       = true;
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=genseg")     == 0)
		 || (strcmp (argv[0], "--format=generalseg") == 0))
			{ lzParams->deGapifyOutput = true;  goto be_general; }

		if ((strcmp (argv[0], "--format=genseg-")     == 0)
		 || (strcmp (argv[0], "--format=generalseg-") == 0))
			{ lzParams->deGapifyOutput = true;  goto be_general_no_header; }

		if ((strcmp (argv[0], "--format=gen")     == 0)
		 || (strcmp (argv[0], "--format=GEN")     == 0)
		 || (strcmp (argv[0], "--format=general") == 0)
		 || (strcmp (argv[0], "--format=GENERAL") == 0))
			{
		be_general:
			free_if_valid ("lzParams->outputInfo", lzParams->outputInfo);
			lzParams->outputFormat = fmtGenpaf;
			lzParams->outputInfo   = copy_string (genpafStandardKeys);
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=gen-")     == 0)
		 || (strcmp (argv[0], "--format=GEN-")     == 0)
		 || (strcmp (argv[0], "--format=general-") == 0)
		 || (strcmp (argv[0], "--format=GENERAL-") == 0))
			{
		be_general_no_header:
			free_if_valid ("lzParams->outputInfo", lzParams->outputInfo);
			lzParams->outputFormat = fmtGenpafNoHeader;
			lzParams->outputInfo   = copy_string (genpafStandardKeys);
			goto next_arg;
			}

		if ((strcmp_prefix (argv[0], "--format=genseg:")     == 0)
		 || (strcmp_prefix (argv[0], "--format=generalseg:") == 0))
			{ lzParams->deGapifyOutput = true;  goto parse_general; }

		if ((strcmp_prefix (argv[0], "--format=genseg-:")     == 0)
		 || (strcmp_prefix (argv[0], "--format=generalseg-:") == 0))
			{ lzParams->deGapifyOutput = true;  goto parse_general_no_header; }

		if ((strcmp_prefix (argv[0], "--format=gen:")     == 0)
		 || (strcmp_prefix (argv[0], "--format=GEN:")     == 0)
		 || (strcmp_prefix (argv[0], "--format=general:") == 0)
		 || (strcmp_prefix (argv[0], "--format=GENERAL:") == 0))
			{
		parse_general:
			free_if_valid ("lzParams->outputInfo", lzParams->outputInfo);
			scan = strchr(argv[0],':') + 1;
			if (*scan == 0)
				suicidef ("empty keys string for --format=general:");
			lzParams->outputFormat = fmtGenpaf;
			lzParams->outputInfo   = parse_genpaf_keys (scan);
			goto next_arg;
			}

		if ((strcmp_prefix (argv[0], "--format=gen-:")     == 0)
		 || (strcmp_prefix (argv[0], "--format=GEN-:")     == 0)
		 || (strcmp_prefix (argv[0], "--format=general-:") == 0)
		 || (strcmp_prefix (argv[0], "--format=GENERAL-:") == 0))
			{
		parse_general_no_header:
			free_if_valid ("lzParams->outputInfo", lzParams->outputInfo);
			scan = strchr(argv[0],':') + 1;
			if (*scan == 0)
				suicidef ("empty keys string for --format=general-:");
			lzParams->outputFormat = fmtGenpafNoHeader;
			lzParams->outputInfo   = parse_genpaf_keys (scan);
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=mapping") == 0)
		 || (strcmp (argv[0], "--format=MAPPING") == 0))
			{
			free_if_valid ("lzParams->outputInfo", lzParams->outputInfo);
			lzParams->outputFormat = fmtGenpaf;
			lzParams->outputInfo   = copy_string (genpafMappingKeys);
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=mapping-") == 0)
		 || (strcmp (argv[0], "--format=MAPPING-") == 0))
			{
			free_if_valid ("lzParams->outputInfo", lzParams->outputInfo);
			lzParams->outputFormat = fmtGenpafNoHeader;
			lzParams->outputInfo   = copy_string (genpafMappingKeys);
			goto next_arg;
			}

		if ((strcmp (argv[0], "--format=rdotplot") == 0))
			{
			free_if_valid ("lzParams->outputInfo", lzParams->outputInfo);
			lzParams->outputFormat   = fmtGenpafNameHeader;
			lzParams->outputInfo     = copy_string (genpafRDotPlotKeys);
			lzParams->deGapifyOutput = true;
			formatIsDotPlot          = true;
			goto next_arg;
			}

		if (strcmp (argv[0], "--format=text") == 0)
			{ lzParams->outputFormat = fmtText;  goto next_arg; }

		if ((strcmp (argv[0], "--format=ztext")    == 0)
		 || (strcmp (argv[0], "--format=zerotext") == 0))
			{ lzParams->outputFormat = fmtZeroText;  goto next_arg; }

		if (strcmp (argv[0], "--format=comp") == 0)
			{ lzParams->outputFormat = fmtHspComp;  goto next_arg; }

		if ((strcmp (argv[0], "--format=diff")        == 0)
		 || (strcmp (argv[0], "--format=diffs")       == 0)
		 || (strcmp (argv[0], "--format=difference")  == 0)
		 || (strcmp (argv[0], "--format=differences") == 0))
			{ lzParams->outputFormat = fmtDiffs;  goto next_arg; }

		if ((strcmp (argv[0], "--format=diff-")        == 0)
		 || (strcmp (argv[0], "--format=diffs-")       == 0)
		 || (strcmp (argv[0], "--format=difference-")  == 0)
		 || (strcmp (argv[0], "--format=differences-") == 0))
			{ lzParams->outputFormat = fmtDiffsNoBlocks;  goto next_arg; }

		if ((strcmp (argv[0], "--format=istats")   == 0)
		 || (strcmp (argv[0], "--format=infstats") == 0))
			{
			if (haveMaxIdentity) goto set_inf_stats_format;
			maxPctId = 70.0;
			goto set_inf_stats_info;
			}

		if (((strcmp_prefix (argv[0], "--format=istats(")   == 0)
		  || (strcmp_prefix (argv[0], "--format=infstats(") == 0))
		 && (argv[0][strlen(argv[0])-1] == ')'))
			{
			scan = strchr(argv[0],'(') + 1;
			scanned = -1;
			sscanf (scan, "%f%n", &maxPctId, &scanned);
			if (scanned == -1) goto cant_understand;
			scan += scanned;
			if (*scan == '%') scan++;
			if (strcmp (scan, ")") != 0) goto cant_understand;
			if ((maxPctId < 0) || (maxPctId > 100)) goto cant_understand;
		set_inf_stats_info:
			lzParams->maxIdentity  = maxPctId / 100.0;
			haveMaxIdentity        = true;
		set_inf_stats_format:
			lzParams->outputFormat = fmtInfStats;
			goto next_arg;
			}

		if (strcmp (argv[0], "--format=identity") == 0)
			{ lzParams->outputFormat = fmtIdDist;  goto next_arg; }

		if (strcmp (argv[0], "--format=deseed") == 0)
			{ lzParams->outputFormat = fmtDeseed;  goto next_arg; }

		if (strcmp (argv[0], "--format=none") == 0)
			{ lzParams->outputFormat = fmtNone;  goto next_arg; }

		if (strcmp (argv[0], "--markend") == 0)
			{ lzParams->endComment = true;  goto next_arg; }

		// --rdotplot=<file>

		if (strcmp_prefix (argv[0], "--rdotplot=") == 0)
			{
			if (lzParams->dotplotFilename != NULL)
				goto duplicated_option;
			lzParams->dotplotFilename = copy_string (argStr);
			goto next_arg;
			}

		// --[no]laj

		if (strcmp (argv[0], "--laj") == 0)
			{ lzParams->lajCompatible = true;  goto next_arg; }

		if (strcmp (argv[0], "--nolaj") == 0)
			{ lzParams->lajCompatible = false;  goto next_arg; }

		// (unadvertised) --expand=<bp>

		if (strcmp_prefix (argv[0], "--expand=") == 0)
			{
			tempInt = string_to_int (argStr);
			if (tempInt < 0)
				suicidef ("--expand cannot be negative");
			else if (tempInt >= 1000)
				suicidef ("--expand must be less than 1000");
			lzParams->textContext = (u32) tempInt;
			goto next_arg;
			}

		// additional file actions

		if ((strcmp_prefix (argv[0], "--action1=[") == 0)
		 && (strcmp_suffix (argv[0], "]")           == 0))
			{
			scan = concatenate_strings (seq1Actions, argStr);
			free_if_valid ("parse_options_loop (seq1Actions)", seq1Actions);
			seq1Actions = scan;
			goto next_arg;
			}

		if ((strcmp_prefix (argv[0], "--action2=[") == 0)
		 && (strcmp_suffix (argv[0], "]")           == 0))
			{
			scan = concatenate_strings (seq2Actions, argStr);
			free_if_valid ("parse_options_loop (seq2Actions)", seq2Actions);
			seq2Actions = scan;
			goto next_arg;
			}

		// precanned expansion arguments

		for (ix=0 ; ix<numExpanders ; ix++)
			{
			if (strcmp (argv[0], expanders[ix].argName) != 0) continue;
			parse_options_string (expanders[ix].expansion, lzParams, izParams);
			goto next_arg;
			}

		// --verbosity=<level>

		if (strcmp_prefix (argv[0], "v=0") == 0)
			{ lzParams->verbosity = 0;  goto next_arg; }

		if (strcmp_prefix (argv[0], "v=1") == 0)
			{ lzParams->verbosity = 10;  goto next_arg; }

		if (strcmp_prefix (argv[0], "--verbosity=") == 0)
			{
			lzParams->verbosity = string_to_int (argStr);
			if      (lzParams->verbosity < 0 ) lzParams->verbosity = 0;
			else if (lzParams->verbosity > 10) lzParams->verbosity = 10;
			goto next_arg;
			}

		if (strcmp_prefix (argv[0], "--gexverbosity=") == 0) // (unadvertised)
			{
			gappedExtendVerbosity = string_to_int (argStr);
			if      (gappedExtendVerbosity < 0 ) gappedExtendVerbosity = 0;
			else if (gappedExtendVerbosity > 10) gappedExtendVerbosity = 10;
			goto next_arg;
			}

		// --[no]runtime

		if (strcmp (argv[0], "--runtime") == 0)
			{ lzParams->reportTiming = true;  goto next_arg; }

		if (strcmp (argv[0], "--noruntime") == 0)
			{ lzParams->reportTiming = false;  goto next_arg; }

		// --tableonly[=count] unadvertised variants

		if (strcmp (argv[0], "--tableonly") == 0)
			{
			lzParams->doSeedSearch = false;
			lzParams->showPosTable = spt_table;
			goto next_arg;
			}

		if (strcmp (argv[0], "--tableonly=count") == 0)
			{
			lzParams->doSeedSearch = false;
			lzParams->showPosTable = spt_countsonly;
			goto next_arg;
			}

		if (strcmp (argv[0], "--tableonly=andcount") == 0)
			{
			lzParams->doSeedSearch = false;
			lzParams->showPosTable = spt_withcounts;
			goto next_arg;
			}

		if (strcmp (argv[0], "--tableonly=distribution") == 0)
			{
			lzParams->doSeedSearch = false;
			lzParams->showPosTable = spt_distribution;
			goto next_arg;
			}

		if (strcmp (argv[0], "--tableonly=stop") == 0)
			{ // (for speed comparisons vs other --tableonly settings)
			lzParams->doSeedSearch = false;
			goto next_arg;
			}

		if (strcmp (argv[0], "--showtable") == 0)
			{ lzParams->showPosTable = spt_table;  goto next_arg; }

		if (strcmp (argv[0], "--showtable=count") == 0)
			{ lzParams->showPosTable = spt_countsonly;  goto next_arg; }

		// --writecapsule=<file>

		if (strcmp_prefix (argv[0], "--writecapsule=") == 0)
			{
			if (lzParams->writeCapsule)
				goto duplicated_option;
			if (lzParams->capsuleFilename != NULL)
				chastise ("can't use --writecapsule with --targetcapsule\n");
			lzParams->capsuleFilename = copy_string (argStr);
			lzParams->writeCapsule    = true;
			lzParams->doSeedSearch    = false;
			goto next_arg;
			}

		// --targetcapsule=<file>

		if (strcmp_prefix (argv[0], "--targetcapsule=") == 0)
			{
			if (lzParams->readCapsule)
				goto duplicated_option;
			if (lzParams->capsuleFilename != NULL)
				chastise ("can't use --targetcapsule with --writecapsule\n");
			if (lzParams->seq1Filename != NULL)
				{
				if (lzParams->seq2Filename != NULL)
					chastise ("can't use --targetcapsule with two queries\n");
				lzParams->seq2Filename = lzParams->seq1Filename;
				lzParams->seq1Filename = NULL;
				}
			lzParams->capsuleFilename = copy_string (argStr);
			lzParams->readCapsule     = true;
			goto next_arg;
			}

		// --[no]stats[=<file>] or (unadvertised) --stats=

		if (strcmp (argv[0], "--stats=") == 0)
			{ lzParams->reportStats = true;  goto next_arg; }

		if (strcmp (argv[0], "--stats") == 0)
			{ lzParams->showStats = true;  goto next_arg; }

		if (strcmp (argv[0], "--nostats") == 0)
			{ lzParams->showStats = false;  goto next_arg; }

		if (strcmp_prefix (argv[0], "--stats=") == 0)
			{
			if (lzParams->statsFilename != NULL) goto duplicated_option;
			lzParams->statsFilename = copy_string (argStr);
			lzParams->showStats     = true;
			goto next_arg;
			}

		// --segments=<file>

		if ((strcmp_prefix (argv[0], "--segments=") == 0)
		 || (strcmp_prefix (argv[0], "--anchors=")  == 0)) // (old option)
			{
			if (lzParams->anchorsFilename != NULL) goto duplicated_option;
			lzParams->anchorsFilename = copy_string (argStr);
			goto next_arg;
			}

		// --notruncationreport (unadvertised)

		if (strcmp (argv[0], "--notruncationreport") == 0)
			{ gapped_extend_inhibitTruncationReport = true;  goto next_arg; }

		// --version

		if ((strcmp (argv[0], "-v")        == 0)
		 || (strcmp (argv[0], "--version") == 0))
			{
			fprintf (helpout, "%s (version %s.%s.%s built %s)\n",
			                  programName,
			                  programVersionMajor, programVersionMinor, programVersionSubMinor, programRevisionDate);

			#if (scoreType == 'I')
			fprintf (helpout, "  score=int");
			#elif (scoreType == 'F')
			fprintf (helpout, " score=float");
			#elif (scoreType == 'D')
			fprintf (helpout, " score=double-float");
			#endif

			fprintf (helpout, ", sequence=%d-bit", maxSequenceIndex);
			fprintf (helpout, ", alloc=%d-bit", maxMallocIndex);

			fprintf (helpout, "\n");
#ifdef __GNUC__
			fprintf (helpout, "  built with gcc-%d.%d.%d\n",
			                  __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif
			exit (EXIT_FAILURE);
			}

		// --help

		if ((strcmp (argv[0], "--all")  == 0)
		 || (strcmp (argv[0], "--help") == 0)
		 || (strcmp (argv[0], "--help=all") == 0))
			{ all_options(); }

		// --help=files

		if (strcmp (argv[0], "--help=files") == 0)
			{ file_options(); }

		// --help=format[s]

		if ((strcmp (argv[0], "--help=format") == 0)
		 || (strcmp (argv[0], "--help=formats") == 0))
			{ format_options(); }

		// --help=short[cuts]

		if ((strcmp (argv[0], "--short") == 0)
		 || (strcmp (argv[0], "--shortcuts") == 0)
		 || (strcmp (argv[0], "--help=short") == 0)
		 || (strcmp (argv[0], "--help=shortcuts") == 0)
		 || (strcmp (argv[0], "--blastz") == 0)
		 || (strcmp (argv[0], "--help=blastz") == 0))
			{ shortcuts(); }

		// --help=yasra

		if ((strcmp (argv[0], "--yasra") == 0)
		 || (strcmp (argv[0], "--help=yasra") == 0))
			{ expander_options ("yasra-specific options", "--yasra"); }

		// --tryout=<what> (unadvertised)

#ifdef tryout
		if (strcmp (argv[0], "--tryout=immediategapped") == 0)
			{
			lzParams->hspToGapped = true;
			goto next_arg;
			}

		if (strcmp_prefix (argv[0], "--tryout=searchlimit:") == 0)
			{
			tempInt = string_to_int (strchr(argv[0],':')+1);
			if (tempInt <= 0)
				suicidef ("--searchlimit must be positive");
			lzParams->searchLimit = tempInt;
			goto next_arg;
			}
#endif // tryout

		// --debug=<what> (unadvertised)

		if (strcmp (argv[0], "--debug") == 0)
			{ debug = 100;  goto next_arg; }

		if (strcmp (argv[0], "--debug=scorematrix") == 0)
			{ dbgShowMatrix = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=sequence") == 0)
			{ sequences_dbgDumpSequence = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=color") == 0)
			{ sequences_dbgAllowColors = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=rawhits") == 0)
			{ seed_search_dbgShowRawHits = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=words") == 0)
			{ pos_table_dbgShowWords = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=maxwordcount") == 0)
			{ pos_table_dbgShowDiscards = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=seedhits") == 0)
			{ seed_search_dbgShowHits = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=seedbases") == 0)
			{ seed_search_dbgShowCoverage = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=chaining") == 0)
			{ chain_dbgChaining = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=chainingtree") == 0)
			{ chain_dbgDumpTree = true;  goto next_arg; }

#ifdef densityFiltering
		if (strcmp (argv[0], "--debug=density") == 0)
			{ seed_search_dbgShowRejections = true;  goto next_arg; }
#endif // densityFiltering

#ifdef snoopHspSubrange
		if (strcmp (argv[0], "--debug=subhsp") == 0)
			{ seed_search_dbgSubrangeHsps = true;  goto next_arg; }
#endif // snoopHspSubrange

		if (strcmp (argv[0], "--debug=currParams") == 0)
			{ dbgShowParams = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=hsps") == 0)
			{ dbgShowHsps = true;  goto next_arg; }

		if ((strcmp (argv[0], "--debug=segments") == 0)
		 || (strcmp (argv[0], "--debug=anchors")  == 0))
			{ dbgShowAnchors = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=sort:diag") == 0)
			{ dbgSortAnchorsByDiag = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=reduction") == 0)
			{ dbgInhibitSegmentReduction = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=masking") == 0)
			{ dbgMasking = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=pctid") == 0)
			{
			gapped_extend_dbgShowIdentity = true;
			identity_dist_dbgShowIdentity = true;
			infer_scores_dbgShowIdentity  = true;
			goto next_arg;
			}

		if (strcmp (argv[0], "--debug=qtobest") == 0)
			{ dna_utilities_dbgShowQToBest = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=qball") == 0)
			{ quantum_dbgQuantumBall = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=maf:diag") == 0)
			{ maf_dbgReportDiag = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=text:diag") == 0)
			{ text_align_dbgReportDiag = true;  goto next_arg; }

		if ((strcmp (argv[0], "--progress") == 0)
		 || (strcmp (argv[0], "--debug=queryprogress") == 0))
			{ dbgQueryProgress = 1;  goto next_arg; }

		if (strcmp_prefix (argv[0], "--progress=") == 0)
			{
			dbgQueryProgress = string_to_int (argStr);
			goto next_arg;
			}

		if (strcmp_prefix (argv[0], "--debug=queryprogress=") == 0)
			{
			dbgQueryProgress = string_to_int (strchr(argStr,'=')+1);
			goto next_arg;
			}

		if ((strcmp (argv[0], "--debug=converge") == 0)
		 || (strcmp (argv[0], "--debug=convergence") == 0))
			{ infer_scores_watchConverge = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=showinferparams") == 0)
			{ infer_scores_showParams = true;  goto next_arg; }

		if (strcmp (argv[0], "--debug=lav+infer") == 0)
			{ infer_scores_outputLav = true;  goto next_arg; }

#ifdef tryout
		if (strcmp (argv[0], "--debug=triviality") == 0)
			{ gapped_extend_dbgTriviality = true;  goto next_arg; }
#endif // tryout

		if (strcmp_prefix (argv[0], "--debug=") == 0)
			{
			debug = string_to_int (argStr);
			if      (debug < 0  ) debug = 0;
			else if (debug > 100) debug = 100;
			goto next_arg;
			}

		// unknown -- argument

		if (strcmp_prefix (argv[0], "--") == 0)
			goto cant_understand;

		// target file name

		if ((lzParams->seq1Filename == NULL) && (!lzParams->readCapsule))
			{
			lzParams->seq1Filename = copy_string (argv[0]);
			goto next_arg_no_copy;
			}

		// query file name

		if (lzParams->seq2Filename == NULL)
			{
			lzParams->seq2Filename = copy_string (argv[0]);
			if (argv[0][0] == '[')
				waywardBracketArg = argv[0];
			goto next_arg_no_copy;
			}

	cant_understand:
		if (strcmp_prefix (argv[0], "--") == 0)
			chastise ("Can't understand \"%s\"\n", argv[0]);
		else if (argv[0][0] == '[')
			chastise ("Can't understand \"%s\"\n"
			          "(my guess) don't use a space between sequence file and bracketed \"action list\"\n",
			          argv[0]);
		else if (waywardBracketArg != NULL)
			chastise ("Can't understand \"%s\"\n"
			          "(my guess) don't use a space between sequence file and %s\n",
			          argv[0], waywardBracketArg);
		else
			chastise ("Can't understand \"%s\"\n"
			          "(my guess) perhaps you have too many sequence files in the command\n",
			          argv[0]);

	duplicated_option:
		chastise ("Duplicated or conflicting option \"%s\"\n", argv[0]);

	next_arg_no_copy:
		if (isTopLevel) lzParams->args[argsLen] = 0; // 'erase' the argument
	next_arg:
		argv++;  argc--;
		continue;

	//next_arg2:
	//	argv+=2;  argc-=2;
	//	continue;
		}
	}

// parse_options_string-- parse options from a string

static void parse_options_string
   (char*		s,
	control*	lzParams,
	control*	izParams)
	{
	int			argC, ix;
	char**		argV, *argS, *scan, *argEnd;

	for (scan=skip_whitespace(s),argC=0 ; *scan!=0 ; argC++)
		{
		scan = skip_darkspace  (scan);
		scan = skip_whitespace (scan);
		}

	argV = malloc_or_die ("parse_options_string (argV)", argC * sizeof(char*));
	argS = copy_string (s);

	for (scan=skip_whitespace(argS),ix=0 ; *scan!=0 ; ix++)
		{
		argV[ix] = scan;
		scan = skip_darkspace  (scan);  argEnd = scan;
		scan = skip_whitespace (scan);
		*argEnd = 0;
		}

	parse_options_loop (argC, argV, lzParams, izParams, /*top level*/ false);

	free_if_valid ("parse_options_string (argV)", argV);
	free_if_valid ("parse_options_string (argS)", argS);
	}

// parse_options-- overall options parsing

static void parse_options
   (int			_argc,
	char**		_argv,
	control*	lzParams,
	control*	izParams)
	{
	int			argc;
	char**		argv;
	int			argsLen, ix;
	score		maxScore;
	int			r, c;
	u8			nuc1, nuc2;
	char*		seq1Filename, *seq2Filename;
	char*		tempS;
	score		myUnitScores[4][4];
	exscoreset*	xss;
	u8			rCh, cCh;

	// skip program name

	argv = _argv+1;  argc = _argc - 1;

	//////////
	// set defaults
	//////////

	*izParams = defaultParams;
	*lzParams = defaultParams;

	seq1Actions           = NULL;
	seq2Actions           = NULL;
	seedString            = NULL;
	scoreFilename         = NULL;
	infControlFilename    = NULL;
	haveXDrop             = false;
	haveYDrop             = false;
	haveStep              = false;
	haveGappedOption      = false;
	haveHspThreshold      = false;
	haveGappedThreshold   = false;
	haveInterpThreshold   = false;
	haveEntropicHsp       = false;
	haveBallScore         = false;
	haveWithTrans         = false;
	haveWithTransForMatch = false;
	haveMaxIdentity       = false;
	useUnitScores         = false;
	gappedExtendVerbosity = -1;
	unitMatch             = 1;
	unitMismatch          = -1;
	haveGapOpen           = false;
	haveGapExtend         = false;
	gapOpenStr            = NULL;
	gapExtendStr          = NULL;
	gapOpen               = 0;
	gapExtend             = 0;
	twinsYes              = defaultTwinsYes;
	minGap                = defaultTwinMinGap;
	maxGap                = defaultTwinMaxGap;
	ballScoreFactor       = -1;	// (indicates we have no factor)
	firstSpecialSub       = NULL;
	formatIsSegments      = false;
	formatIsDotPlot       = false;

	// create a string to copy arguments into (this will be slightly bigger
	// than we need)

	argsLen = 2;
	for (ix=0 ; ix<argc ; ix++)
		argsLen += 1 + strlen(argv[ix]);

	lzParams->args = malloc_or_die ("parse_options (lzParams->args)", argsLen);
	lzParams->args[0] = 0;

	for (r=0 ; r<4 ; r++)
			for (c=0 ; c<4 ; c++)
		specialSubScores[r][c] = worstPossibleScore;

	//////////
	// scan arguments
	//////////

	parse_options_loop (argc, argv, lzParams, izParams, /*top level*/ true);

	if (lzParams->outputFilename == NULL)
		lzParams->outputFile = stdout;
	else
		lzParams->outputFile = fopen_or_die (lzParams->outputFilename, "wt");

	if ((lzParams->dotplotFilename) && (formatIsDotPlot))
		suicidef ("--format=rdotplot can't be used with --rdotplot=<file>");

	if (lzParams->dotplotFilename != NULL)
		lzParams->dotplotFile = fopen_or_die (lzParams->dotplotFilename, "wt");

	//////////
	// do some post-processing of sequence names
	//////////

	// bind 'extra' actions to sequence names

	if (seq1Actions != NULL)
		{
		// $$$ this is too stringent-- certain actions would be ok, so we
		// $$$ .. really ought to just check for the actions that would be bad;
		// $$$ .. however, we'd have no target filename string to bind the
		// $$$ .. actions to
		if (lzParams->readCapsule)
			suicidef ("--action1 can't be used with --targetcapsule");
		tempS = concatenate_strings (lzParams->seq1Filename, seq1Actions);
		free_if_valid ("parse_options (seq1Filename)", lzParams->seq1Filename);
		free_if_valid ("parse_options (seq1Actions)",  seq1Actions);  seq1Actions = NULL;
		lzParams->seq1Filename = tempS;
		}

	if (seq2Actions != NULL)
		{
		if (lzParams->seq2Filename == NULL)
			suicidef ("--action2 can't be used without query sequence file");
		tempS = concatenate_strings (lzParams->seq2Filename, seq2Actions);
		free_if_valid ("parse_options (seq2Filename)", lzParams->seq2Filename);
		free_if_valid ("parse_options (seq2Actions)",  seq2Actions);  seq2Actions = NULL;
		lzParams->seq2Filename = tempS;
		}

	// determine which sequences are quantum

	lzParams->targetIsQuantum = false;
	if ((lzParams->seq1Filename != NULL) && (!lzParams->readCapsule))
		lzParams->targetIsQuantum = name_spec_is_quantum (lzParams->seq1Filename);

	if (lzParams->seq2Filename != NULL)
		lzParams->queryIsQuantum = name_spec_is_quantum (lzParams->seq2Filename);

	if ((lzParams->targetIsQuantum) || (lzParams->queryIsQuantum))
		{
		if (lzParams->inferScores)
			suicide ("scoring inference cannot be performed with quantum DNA");
		if ((lzParams->minIdentity > 0) || (lzParams->maxIdentity < 1))
			suicide ("identity filtering cannot be used with quantum DNA");
		if (lzParams->minMatchCount > 0)
			suicide ("match count filtering cannot be used with quantum DNA");
		if (lzParams->outputFormat == fmtIdDist)
			suicide ("--format=identity cannot be used with quantum DNA");
#ifdef densityFiltering
		if (lzParams->maxDensity != 0)
			suicide ("--density cannot be used with quantum DNA");
#endif // densityFiltering
		}

	//////////
	// check for sensibility
	//////////

#ifdef disallowEntropy
	if (lzParams->entropicHsp)
		chastise ("--entropy is currently disabled\n");
#endif

#ifndef collect_stats
	if (lzParams->showStats)
		chastise ("--stats is not implemented in this build of the program\n");
#endif // collect_stats

	if (lzParams->writeCapsule)
		lzParams->outputFormat = fmtNone;

	if ((lzParams->seq1Filename == NULL) && (!lzParams->readCapsule))
		chastise ("You must specify a target file\n");

	if (lzParams->inferOnly)
		{
		if (lzParams->noHitFiltering)
			chastise ("--rawhits can't be used with --inferonly\n");

		if (lzParams->dynamicMasking > 0)
			chastise ("--masking can't be used with --inferonly\n");

		if (lzParams->reportCensus)
			chastise ("--census can't be used with --inferonly\n");

		if (lzParams->outputFormat != defaultParams.outputFormat)
			chastise ("--format=%s can't be used with --inferonly\n",
			          formatNames[lzParams->outputFormat]);

		if (lzParams->innerThreshold > 0)
			chastise ("--inner can't be used with --inferonly\n");

		if (lzParams->anchorsFilename != NULL)
			chastise ("--segments can't be used with --inferonly\n");
		}

	if (lzParams->selfCompare)
		{
		if (lzParams->seq2Filename != NULL)
			chastise ("--self can't be used when you specify a query file\n");
		if (lzParams->anchorsFilename != NULL)
			chastise ("--segments can't be used with --self\n");
		lzParams->seq2Filename = copy_string (lzParams->seq1Filename);
		if (lzParams->mirror == -1)
			lzParams->mirror = true;
		}
	else if (lzParams->mirror == true)
		chastise ("--mirror can only be used with --self\n");
	else
		lzParams->mirror = false;

	if (lzParams->readCapsule)
		{
		capseed* seedCapsule;

		if (seedString != NULL)
			chastise ("can't set word size or seed pattern with --targetcapsule\n");
		if (haveStep)
			chastise ("can't use --step with --targetcapsule\n");
		if (lzParams->dynamicMasking > 0)
			chastise ("can't use --masking with --targetcapsule\n");
		if (lzParams->wordCountLimit > 0)
			chastise ("can't use --maxwordcount with --targetcapsule\n");
		if (lzParams->maxIndexBits != defaultParams.maxIndexBits)
			chastise ("can't use --word with --targetcapsule\n");

		lzParams->capsule = open_capsule_file (lzParams->capsuleFilename);
		seedCapsule = locate_capsule_data (lzParams->capsule, cap_seed,
		                                   NULL, NULL);
		if (seedCapsule == NULL)
			suicide ("bad capsule file (missing seed)");
		lzParams->step = seedCapsule->step;
		}

	if (lzParams->writeCapsule)
		{
		if (lzParams->seq2Filename != NULL)
			chastise ("--writecapsule can't be used when you specify a query file\n");
		if (lzParams->inferScores)
			chastise ("can't use --infer with --writecapsule\n");
		if (lzParams->anchorsFilename != NULL)
			chastise ("can't use --segments with --writecapsule\n");
		if (haveXDrop)
			chastise ("can't use --xdrop with --writecapsule\n");
		if (haveYDrop)
			chastise ("can't use --ydrop with --writecapsule\n");
		if (haveHspThreshold)
			chastise ("can't use --hspthresh with --writecapsule\n");
		if (haveGappedThreshold)
			chastise ("can't use --gappedthresh with --writecapsule\n");
		if (haveInterpThreshold)
			chastise ("can't use --inner with --writecapsule\n");
		if (haveEntropicHsp)
			chastise ("can't use --entropy with --writecapsule\n");
		if (haveBallScore)
			chastise ("can't use --ball with --writecapsule\n");
		if ((haveWithTrans) && (!haveWithTransForMatch))
			chastise ("can't use --transition with --writecapsule\n");
		if (haveMaxIdentity)
			chastise ("can't use --identity with --writecapsule\n");
		if ((haveGapOpen) || (haveGapExtend))
			chastise ("can't use --gap with --writecapsule\n");
		}

	if (!lzParams->doSeedSearch)
		{
		if (lzParams->seq2Filename != NULL)
			chastise ("--tableonly can't be used when you specify a query file\n");
		if (lzParams->inferScores)
			chastise ("--infer and --tableonly are not compatible\n");
		}

	if (lzParams->maxIndexBits < 8)
		chastise ("--word doesn't allow so few bits (%d)\n",
				  lzParams->maxIndexBits);

	if (lzParams->tracebackMem < 100*1024)
		chastise ("--traceback must be at least 100K (it's only %s)\n",
				  unitize(lzParams->tracebackMem));

#ifndef noSeedHitQueue
	if (lzParams->seedHitQueueSize < 0)
		chastise ("--seedqueue can't be negative\n");
#endif // not noSeedHitQueue

	if ((lzParams->reportCensus)
	 && (lzParams->censusFilename == NULL))
		{
		if ((lzParams->outputFormat != fmtLav)
		 && (lzParams->outputFormat != fmtLavComment)
		 && (lzParams->outputFormat != fmtLavScore)
		 && (lzParams->outputFormat != fmtLavText))
			chastise ("--census with --format=%s requires --census=<file>\n",
			          formatNames[lzParams->outputFormat]);
		}

	if (lzParams->hspToGapped)
		{
		// $$$ hspToGapped should turn off second-stage gapped, although nothing
		//     .. will happen in the second stage, since the anchors list will
		//     .. be empty

		if (lzParams->inferScores)
			chastise ("can't use --anyornone with --infer[only]\n");

		if (lzParams->innerThreshold > 0)
			chastise ("can't use --anyornone with --inner\n");

		if (lzParams->anchorsFilename != NULL)
			chastise ("can't use --anyornone with --segments\n");

		if (lzParams->hspThreshold.t != 'S')
			chastise ("can't use --anyornone with adaptive hsp score threshold\n");

		if (lzParams->chain)
			chastise ("can't use --anyornone with --chain\n");

		if (lzParams->innerThreshold > 0)
			chastise ("can't use --anyornone with --inner\n");
		}

	if (lzParams->searchLimit > 0)
		{
		if (lzParams->inferScores)
			chastise ("can't use --anyornone with --infer[only]\n");

		if (lzParams->innerThreshold > 0)
			chastise ("can't use --anyornone with --inner\n");

		if (lzParams->anchorsFilename != NULL)
			chastise ("can't use --anyornone with --segments\n");

		if (lzParams->hspThreshold.t != 'S')
			chastise ("can't use --anyornone with adaptive hsp score threshold\n");

		if ((lzParams->targetIsQuantum) || (lzParams->queryIsQuantum))
			chastise ("can't use --anyornone with quantum dna\n");
		}

	if ((formatIsSegments) && (!haveGappedOption))
		{
		if (haveInterpThreshold)
			chastise ("--inner cannot be used with --writesegments\n");
		lzParams->gappedExtend = false;
		}

	//////////
	// set up score set
	//////////

	if (lzParams->inferScores)
		{
#if ((scoreType == 'I') && (!defined infer_anything))
		suicide ("scoring inference can't be performed with integer arithmetic;  use lastz_D");
#endif

		if (lzParams->anchorsFilename != NULL)
			chastise ("--segments can't be used with --infer[only]\n");

		if (scoreFilename != NULL)
			chastise ("can't use --infer[only] and --scores together\n");

		if (useUnitScores)
			chastise ("can't use --infer[only] and --match (or --unitscores) together\n");

		if (haveGapOpen)
			chastise ("can't use --infer[only] and --gap (or O=) together\n");

		if (haveGapExtend)
			chastise ("can't use --infer[only] and --gap (or E=) together\n");

		if (firstSpecialSub != NULL)
			chastise ("can't use --infer[only] and special substitution scores together\n");
		}

	if (lzParams->gfExtend == gfexNoExtend)
		{
		if ((!lzParams->gappedExtend)
		 && (scoreFilename != NULL)
		 && (!lzParams->targetIsQuantum)
		 && (!lzParams->queryIsQuantum))
			chastise ("--scores requires --gfextend or --gapped\n");

		if (haveXDrop)
			chastise ("--xdrop requires --gfextend\n");

		if (haveHspThreshold)
			chastise ("--hspthresh requires --gfextend\n");

		if (haveEntropicHsp)
			chastise ("--entropy requires --gfextend\n");

		if (lzParams->xDropUntrimmed)
			chastise ("--noxtrim requires --gfextend\n");

		lzParams->xDrop          = 0;
		lzParams->hspThreshold.t = 'S';
		lzParams->hspThreshold.s = 0;
		lzParams->entropicHsp    = false;
		}

	if (!lzParams->chain)
		{
		if (lzParams->chainDiag != 0)
			chastise ("G=<score> requires --chain\n");

		if (lzParams->chainAnti != 0)
			chastise ("R=<score> requires --chain\n");
		}

	if (lzParams->chain)
		{
		if (lzParams->anchorsFilename != NULL)
			chastise ("--segments can't be used with --chain\n");
		}

	if (!lzParams->gappedExtend)
		{
		if ((haveGapOpen) || (haveGapExtend))
			chastise ("--gap (or O= or E=) requires --gapped\n");

		if (haveYDrop)
			chastise ("--ydrop requires --gapped\n");

		if (haveGappedThreshold)
			chastise ("--gappedThreshold requires --gapped\n");

		if (haveInterpThreshold)
			chastise ("--inner requires --gapped\n");

		if (lzParams->yDropUntrimmed)
			chastise ("--noytrim requires --gapped\n");

		if ((lzParams->maxContinuity < 1)
		 && (!lzParams->doSeedSearch)
		 && (!lzParams->writeCapsule))
			chastise ("--continuity maximum less than 1 requires --gapped\n");

		if (lzParams->gappedAllBounds)
			chastise ("--allgappedbounds requires --gapped\n");
		}

	if (lzParams->gappedExtend)
		{
		if (formatIsSegments)
			chastise ("can't used --writesegments with --gapped\n");
		}

	if (lzParams->anchorsFilename != NULL)
		{
		if (haveHspThreshold)
			chastise ("--segments can't be used with --hspthresh\n");
		if (haveXDrop)
			chastise ("--segments can't be used with --xdrop\n");
		if (seedString != NULL)
			chastise ("can't set word size or seed pattern with --segments\n");
		}

	if (haveXDrop && (lzParams->xDrop <= 0))
		chastise ("%d is not a legitimate x-drop threshold\n", lzParams->xDrop);

	if (haveYDrop && (lzParams->yDrop <= 0))
		chastise ("%d is not a legitimate y-drop threshold\n", lzParams->yDrop);

	if ((useUnitScores) && (scoreFilename != NULL))
		chastise ("can't use --match (or --unitscores) and --scores together\n");

	if (scoreFilename != NULL)
		{
		// read scores and score-related paramters from file, allowing them
		// to be overridden by the command line

		xss = read_score_set_by_name (scoreFilename);
		lzParams->scoring = (scoreset*) xss;

		if (xss->seedSet)
			{
			if (seedString == NULL) // it contains params in command-line syntax
				parse_options_loop (1, &xss->seed, lzParams, izParams, /*top level*/ false);
			free_if_valid ("xss->seed", xss->seed);  xss->seed = NULL;
			}

		if (!haveGapOpen)
			gapOpen = lzParams->scoring->gapOpen;
		else
			{
			lzParams->scoring->gapOpen    = gapOpen;
			lzParams->scoring->gapOpenSet = true;
			}

		if (!haveGapExtend)
			gapExtend = lzParams->scoring->gapExtend;
		else
			{
			lzParams->scoring->gapExtend    = gapExtend;
			lzParams->scoring->gapExtendSet = true;
			}

		if ((!haveHspThreshold) && (xss->hspThresholdSet))
			{
			lzParams->hspThreshold.t    = 'S';
			lzParams->hspThreshold.s    = xss->hspThreshold;
			haveHspThreshold            = true;
			}
		if ((!haveGappedThreshold) && (xss->gappedThresholdSet))
			{
			lzParams->gappedThreshold.t = 'S';
			lzParams->gappedThreshold.s = xss->gappedThreshold;
			haveGappedThreshold         = true;
			}
		if ((!haveXDrop) && (xss->xDropSet))
			{
			lzParams->xDrop             = xss->xDrop;
			haveXDrop                   = true;
			}
		if ((!haveYDrop) && (xss->yDropSet))
			{
			lzParams->yDrop             = xss->yDrop;
			haveYDrop                   = true;
			}
		if ((!haveBallScore) && (ballScoreFactor < 0) && (xss->ballScoreSet))
			{
			if (xss->ballScoreFactor < 0)
				{ lzParams->ballScore = xss->ballScore;  haveBallScore = true; }
			else
				ballScoreFactor = xss->ballScoreFactor;
			}
		if ((!haveStep) && (xss->stepSet))
			{
			lzParams->step              = xss->step;
			haveStep                    = true;
			}

		if ((haveGapOpen) && (gapOpen + gapExtend < 0))
			chastise ("%s is not a valid gap open penalty\n", gapOpenStr);
		if ((haveGapExtend) && (gapExtend < 0))
			chastise ("%s is not a valid gap extension penalty\n", gapExtendStr);
		}
	else if (useUnitScores)
		{
		// use unit scoring matrix, scaled if requested

		scratchThreshold.t = 'S';
		if (scoreType == 'I')
			scratchThreshold.s = (score) ceil (unitScores_thresh * unitMatch);
		else
			scratchThreshold.s = (score)      (unitScores_thresh * unitMatch);

		if (!haveGapOpen)
			{
			if (scoreType == 'I')
				gapOpen = (score) ceil (unitScores_open * -unitMismatch);
			else
				gapOpen = (score)      (unitScores_open * -unitMismatch);
			haveGapOpen = true;
			}
		if (!haveGapExtend)
			{
			if (scoreType == 'I')
				gapExtend = (score) ceil (unitScores_extend * -unitMismatch);
			else
				gapExtend = (score)      (unitScores_extend * -unitMismatch);
			haveGapExtend = true;
			}

		if (!haveHspThreshold)
			{
			lzParams->hspThreshold = scratchThreshold;
			haveHspThreshold = true;
			}

		if ((!haveGappedThreshold) && (lzParams->gfExtend == gfexExact))
			{
			lzParams->gappedThreshold = scratchThreshold;
			haveGappedThreshold = true;
			}

		if ((!haveXDrop) && (!lzParams->inferScores))
			{
			if (scoreType == 'I')
				lzParams->xDrop = (score) ceil (10.0 * sqrt(-unitMismatch));
			else
				lzParams->xDrop = (score)      (10.0 * sqrt(-unitMismatch));
			haveXDrop = true;
			}

		if ((!haveYDrop) && (!lzParams->inferScores))
			{
			lzParams->yDrop = 2 * lzParams->xDrop;
			haveYDrop = true;
			}

		if ((haveGapOpen) && (gapOpen + gapExtend < 0))
			chastise ("%s is not a valid gap open penalty\n", gapOpenStr);
		if ((haveGapExtend) && (gapExtend < 0))
			chastise ("%s is not a valid gap extension penalty\n", gapExtendStr);

		for (r=0 ; r<4 ; r++)
				for (c=0 ; c<4 ; c++)
			myUnitScores[r][c] = (r==c)? unitMatch : unitMismatch;

		lzParams->scoring = new_dna_score_set (myUnitScores,
		                                       unitScores_X    * -unitMismatch,
		                                       unitScores_fill * -unitMismatch,
		                                       gapOpen, gapExtend);
		}
	else if (lzParams->inferScores)
		{
		; // (do nothing, lzParams->scoring will be created by inference)
		}
	else
		{
		// use blastz default scoring matrix
		if (!haveGapOpen)   gapOpen   = HOXD70_open;
		if (!haveGapExtend) gapExtend = HOXD70_extend;
		if ((haveGapOpen) && (gapOpen + gapExtend < 0))
			chastise ("%s is not a valid gap open penalty\n", gapOpenStr);
		if ((haveGapExtend) && (gapExtend < 0))
			chastise ("%s is not a valid gap extension penalty\n", gapExtendStr);
		lzParams->scoring = new_dna_score_set (HOXD70,
		                                      HOXD70_X, HOXD70_fill,
		                                      gapOpen, gapExtend);
		}

	if (firstSpecialSub != NULL)
		{
		score worstScore = 0;

		if ((!lzParams->scoring->rowsAreDna)
		 || (!lzParams->scoring->colsAreDna))
			suicidef ("special substitution scores (e.g. %s) can't be used with quantum DNA scores",
			          firstSpecialSub);

		free_if_valid ("parse_options (firstSpecialSub)", firstSpecialSub);

		for (r=0 ; r<4 ; r++)
				for (c=0 ; c<4 ; c++)
			{
			if (specialSubScores[r][c] == worstPossibleScore)
				{
				nuc1 = bits_to_nuc[r];
				nuc2 = bits_to_nuc[c];
				specialSubScores[r][c] = lzParams->scoring->sub[nuc1][nuc2];
				}
			if (specialSubScores[r][c] < worstScore)
				worstScore = specialSubScores[r][c];
			}
		free_score_set ("lzParams->scoring (special subs)", lzParams->scoring);
		lzParams->scoring = new_dna_score_set (specialSubScores,
		                                      10*worstScore, worstScore,
		                                      gapOpen, gapExtend);
		}

	//////////
	// convert seed string to a seed structure
	//////////

	if ((lzParams->targetIsQuantum) || (lzParams->queryIsQuantum))
		{
		if ((haveWithTrans) && (lzParams->withTrans != 0))
			suicidef ("can't use --transitions with quantum DNA",
		              lzParams->seq2Filename);
		lzParams->withTrans = 0;
		}

	create_seed_structure (lzParams, &seedString,
	                       haveWithTrans, twinsYes, minGap, maxGap);

	if ((lzParams->targetIsQuantum) || (lzParams->queryIsQuantum))
		{
		if (lzParams->hitSeed->type != 'S')
			suicide ("quantum DNA requires a strict seed\n"
			          "(only 1s and 0s allowed, no Ts, no --seed=half)");
		}

	if (pos_table_dbgShowDiscards)
		pos_table_dbgSeed = lzParams->hitSeed;

	//////////
	// compute default values for parameters that have not been set by the user
	//////////

	if (!haveXDrop)
		{
		if (lzParams->inferScores)
			lzParams->xDrop = -1;	// (will fill in after scoring inference)
		else
			{
			rCh = lzParams->scoring->rowChars[0];
			cCh = lzParams->scoring->colChars[0];
			lzParams->xDrop = 10 * lzParams->scoring->sub[rCh][cCh];
			}
		}

	if (!haveYDrop)
		{
		if (lzParams->inferScores)
			lzParams->yDrop = -1;	// (will fill in after scoring inference)
		else
			lzParams->yDrop = lzParams->scoring->gapOpen + 300 * lzParams->scoring->gapExtend;
		}

	if (!haveGappedThreshold)
		{
		if (lzParams->gfExtend == gfexExact)
			lzParams->gappedThreshold = defaultParams.hspThreshold;
		else
			lzParams->gappedThreshold = lzParams->hspThreshold;
		}

	//////////
	// set up others
	//////////

//ии somewhere I need to check for a mismatch between strandedness and the
//ии complementabilty of the score set (for quantum query)

	// create a version of the scoring set that penalizes lowercase bases

	if (!lzParams->inferScores)
		lzParams->maskedScoring = masked_score_set (lzParams->scoring);

	// make scores vs N be ambiguous, if desired and if rows and columns are
	// both DNA

	if ((lzParams->nIsAmbiguous)
	 && (!lzParams->scoring->rowsAreDna)
	 && (!lzParams->scoring->colsAreDna))
		suicidef ("can't use --ambiguous if both target and query are quantum");

	if (lzParams->allowAmbiDNA)
		{
		ambiguate_iupac (lzParams->scoring,       0, -2*lzParams->scoring->gapExtend);
		ambiguate_iupac (lzParams->maskedScoring, 0, -2*lzParams->scoring->gapExtend);
		}
	if (lzParams->nIsAmbiguous)
		{
		ambiguate_n (lzParams->scoring,       0, -2*lzParams->scoring->gapExtend);
		ambiguate_n (lzParams->maskedScoring, 0, -2*lzParams->scoring->gapExtend);
		}

	if (dbgShowMatrix)
		{
		fprintf        (stderr, "lzParams->scoring:\n");
		dump_score_set (stderr, lzParams->scoring,       (u8*)"ACGTacgtNnXF", (u8*)"ACGTacgtNnXF");
		fprintf        (stderr, "\n");
		fprintf        (stderr, "lzParams->maskedScoring:\n");
		dump_score_set (stderr, lzParams->maskedScoring, (u8*)"ACGTacgtNnXF", (u8*)"ACGTacgtNnXF");
		}

	if (lzParams->inferScores)
		maxScore = 0;  // (maxScore is not needed)
	else
		maxScore = max_in_score_matrix (lzParams->scoring);

	if (lzParams->seq1Filename != NULL) seq1Filename = lzParams->seq1Filename;
	                               else seq1Filename = "(unnamed target file)";

	if (lzParams->seq2Filename != NULL) seq2Filename = lzParams->seq2Filename;
	                               else seq2Filename = "(unnamed query file)";

	if (!lzParams->inferScores)
		{
		if ((!lzParams->targetIsQuantum)				// DNA target
		 && (!lzParams->maskedScoring->rowsAreDna))
			suicidef ("row scores are for quantum DNA, but target is not");

		if ((lzParams->doSeedSearch)
		 && (!lzParams->queryIsQuantum)					// DNA query
		 && (!lzParams->maskedScoring->colsAreDna))
			suicidef ("column scores are for quantum DNA, but query is not");

		if ((lzParams->targetIsQuantum)					// quantum DNA target
		 && (lzParams->maskedScoring->rowsAreDna))
			suicidef ("target is quantum DNA, but row scores are not");

		if ((lzParams->queryIsQuantum)					// quantum DNA query
		 && (lzParams->maskedScoring->colsAreDna))
			suicidef ("query is quantum DNA, but column scores are not");

		if (((haveBallScore) || (ballScoreFactor >= 0))
		 && ((!lzParams->targetIsQuantum) && (!lzParams->queryIsQuantum)))
			suicidef ("--ball can't be used with DNA target and query");

		if ((lzParams->targetIsQuantum)
		 || (lzParams->queryIsQuantum))
			{
			if ((haveBallScore) && (lzParams->ballScore < 0))
				chastise (scoreFmtSimple " is not a legitimate ball threshold\n",
						  lzParams->ballScore);
			if (!haveBallScore)
				{
				if (ballScoreFactor < 0) ballScoreFactor = defaultBallScoreFactor;
				lzParams->ballScore = ballScoreFactor * maxScore
								   * lzParams->hitSeed->weight/2;
				}
			if (lzParams->ballScore < 0)
				suicidef ("quantum DNA (%s) requires --ball", seq2Filename);
			if ((lzParams->outputFormat == fmtAxt)
			 || (lzParams->outputFormat == fmtAxtComment))
				suicidef ("--axt doesn't support quantum DNA");
			if ((lzParams->outputFormat == fmtMaf)
			 || (lzParams->outputFormat == fmtMafComment)
			 || (lzParams->outputFormat == fmtMafNoComment))
				suicidef ("--maf doesn't support quantum DNA");
			if ((lzParams->outputFormat == fmtGenpaf)
			 && (strchr (lzParams->outputInfo, genpafText1) != NULL))
				suicidef ("--format=general:text1 doesn't support quantum DNA");
			if ((lzParams->outputFormat == fmtGenpaf)
			 && (strchr (lzParams->outputInfo, genpafText2) != NULL))
				suicidef ("--format=general:text2 doesn't support quantum DNA");
			}
		}

	// build a seed for interpolation (an n-mer exact match)

	if (lzParams->innerThreshold > 0)
		{
		char seedString[innerWordSize+1];
		int i;

		for (i=0 ; i<innerWordSize ; i++)
			seedString[i] = '1';
		seedString[innerWordSize] = 0;

		parse_seeds_string (seedString, &lzParams->innerSeed, 28);
		lzParams->innerSeed->withTrans = 0;
		}

	// decide whether we need to 'waste' time reading files just to get the
	// true sequence length (some output formats reuire the sequence length,
	// others don't)

	lzParams->needTrueLengths
	   =  (lzParams->anchorsFilename != NULL)
	   || (lzParams->outputFormat == fmtAxt)
	   || (lzParams->outputFormat == fmtAxtComment)
	   || (lzParams->outputFormat == fmtMaf)
	   || (lzParams->outputFormat == fmtMafComment)
	   || (lzParams->outputFormat == fmtMafNoComment)
	   || (lzParams->outputFormat == fmtSoftSam)
	   || (lzParams->outputFormat == fmtHardSam)
	   || (lzParams->outputFormat == fmtCigar)
	   || (((lzParams->outputFormat == fmtGenpaf)
	     || (lzParams->outputFormat == fmtGenpafNameHeader))
	    && ((strchr (lzParams->outputInfo, genpafSize1) != NULL)
	     || (strchr (lzParams->outputInfo, genpafSize2) != NULL)));

	// propagate some deep control to other modules

	if (lzParams->verbosity >= 2)
		{
		showProgress             = true;
		pos_table_showProgress   = true;
		seed_search_showProgress = true;
		}

	gapped_extend_verbosity      = (gappedExtendVerbosity>=0)? gappedExtendVerbosity
	                                                         : lzParams->verbosity;
	gapped_extend_dbgShowHsps    = dbgShowHsps;
	gapped_extend_dbgShowAnchors = dbgShowAnchors;

	sequences_keepFastaArrow     = lzParams->lajCompatible;

	//////////
	// set up inference parameters
	//
	// note 1: Ideally we want the default gapOpen and gapExtend penalties to
	//         be relative to the worst substitution score.  Unfortunately,
	//         if scores are integers we don't have any mechanism to pass the
	//         ratios in the score variable (the ratios are not integers).  So
	//         the code below does what we want when scores are non-integers,
	//         but makes a poor approximation attempt otherwise.  This could be
	//         corrected with some effort, but perfroming inference with integer
	//         scores is undesirable, for other reason, so the effort is not
	//         warrented.
	//////////

	if (lzParams->inferScores)
		{
		infcontrol tempIc;

		tempIc = izParams->ic;
		*izParams = *lzParams;
		izParams->ic = tempIc;
		izParams->tracebackMem = 0;
		izParams->outputFormat = fmtNone;

		izParams->hitSeed = copy_seeds (lzParams->hitSeed);

		if (izParams->scoring == NULL)
			izParams->scoring = new_dna_score_set (unitScores,
			                                       unitScores_X, unitScores_fill,
			                                       worstPossibleScore, worstPossibleScore);

		if (infControlFilename != NULL)
			read_control_file_by_name (infControlFilename, izParams);

		if ((izParams->ic.inferScale >  0)
		 && (izParams->ic.inferScale != 1))
			scale_score_set (izParams->scoring, izParams->ic.inferScale);

		if (izParams->scoring->gapOpen == worstPossibleScore)
			{ // (see note 1 about 20 lines above)
			if (scoreType != 'I')
				{
				izParams->ic.gapOpenIsRatio = ratioMinSubScore;
				izParams->scoring->gapOpen  = unitScores_open;
				}
			else if (izParams->ic.inferScale > 0)
				izParams->scoring->gapOpen = (score) ceil (unitScores_open * izParams->ic.inferScale);
			else
				izParams->scoring->gapOpen =               unitScores_open;
			}

		if (izParams->scoring->gapExtend == worstPossibleScore)
			{ // (see note 1 about 30 lines above)
			if (scoreType != 'I')
				{
				izParams->ic.gapExtendIsRatio = ratioMinSubScore;
				izParams->scoring->gapExtend  = unitScores_extend;
				}
			else if (izParams->ic.inferScale > 0)
				izParams->scoring->gapExtend = (score) ceil (unitScores_extend * izParams->ic.inferScale);
			else
				izParams->scoring->gapExtend =               unitScores_extend;
			}

		izParams->maskedScoring = masked_score_set (izParams->scoring);
		}

	//////////
	// clean up
	//////////

	free_if_valid ("parse_options (seedString)",         seedString);
	free_if_valid ("parse_options (scoreFilename)",      scoreFilename);
	free_if_valid ("parse_options (infControlFilename)", infControlFilename);
	}

//----------
//
// create_seed_structure--
//	Convert the user-specified parameters into a seed and set any related
//	variables.
//
//----------
//
// Arguments:
//	control*	lzParams:	Control data for the primary alignment.  This
//							.. routine will potnetially modify several fields.
//	char**		seedString:	The user-defined seed pattern.  This may be
//							.. modified by this routine.  Specifically, if the
//							.. string is NULL, a nely allocated string is
//							.. created and copied to this.
//
// Returns:
//	(nothing)
//
//----------

static void create_seed_structure
   (control*	lzParams,
	char**		_seedString,
	int			haveWithTrans,
	int			twinsYes,
	int			minGap,
	int			maxGap)
	{
	char*		seedString = *_seedString;

	// reconstruct seed from the capsule

	if (lzParams->capsule != NULL)
		{
		capseed*	seedCapsule;
		u64			seedCapsuleSize, expectedSize;
		int			shift[100];		// (100 is a safe upper bound;  usually
		u32*		mask;			//  .. fewer than five are needed)
		u32*		transFlips;
		int			numParts, numFlips, ix;
		u32*		scan;

		seedCapsule = locate_capsule_data (lzParams->capsule, cap_seed,
		                                   NULL, &seedCapsuleSize);
		if (seedCapsule == NULL)
			suicide ("bad capsule file (missing seed)");

		numParts = (int) seedCapsule->numParts;
		if (numParts > (int) (sizeof(shift)/sizeof(shift[0])))
			suicidef ("internal error handling capsule file (numParts = %d)",
			          numParts);

		scan       = &seedCapsule->shift0;
		mask       = &scan[numParts];
		transFlips = &mask[numParts];

		for (ix=0 ; ix<numParts ; ix++)
			shift[ix] = scan[ix];

		for (numFlips=0 ; transFlips[numFlips]!=0 ; numFlips++)
			;

		expectedSize = (sizeof(capseed) - sizeof(u32))	// standard fields
		             + (numParts * sizeof(u32))			// shift[] array
		             + (numParts * sizeof(u32))			// mask[] array
		             + ((numFlips+1) * sizeof(u32));	// transFlips[] array
		if (seedCapsuleSize != expectedSize)
			suicidef ("bad capsule file seed (expected 0x%s, actual 0x%s)",
			          hex_64_string(expectedSize), hex_64_string(seedCapsuleSize));

		lzParams->hitSeed = reconstruct_seed
		                      ((char) seedCapsule->type,
		                       (int)  seedCapsule->length,
		                       (int)  seedCapsule->weight,
		                              NULL, // $$$ add support for capsules containing seed pattern
		                              seedCapsule->resolvingMask,
		                       (int)  seedCapsule->revComp,
		                       (int)  seedCapsule->isHalfweight,
		                       (int)  seedCapsule->numParts,
		                              shift, mask, transFlips);
		}

	// create seed from user params

	else
		{
		if (seedString == NULL)
			*_seedString = seedString = copy_string (defaultSeedString);

		parse_seeds_string (seedString, &lzParams->hitSeed, lzParams->maxIndexBits);
		}

	if ((lzParams->hitSeed->type == 'H')
	  && (!haveWithTrans))
		lzParams->withTrans = 0;

	if ((lzParams->withTrans != 0)
	 && (lzParams->hitSeed->type != 'S')
	 && (lzParams->hitSeed->transFlips == NULL))
		chastise ("--transition can only be used with strict seeds (1s and 0s)\n");

	if (lzParams->minMatches >= 0)
		{
		int   numPositions;
		char* pScan;

		if ((lzParams->filterCaresOnly) && (lzParams->hitSeed->pattern == NULL))
			chastise ("--filter=cares: can not be used with a patternless seed\n");

		if (!lzParams->filterCaresOnly)
			numPositions = lzParams->hitSeed->length;
		else
			{
			numPositions = 0;
			for (pScan=lzParams->hitSeed->pattern ; *pScan!=0 ; )
				{ if (*(pScan++) != '0') numPositions++; }
			}

		if (lzParams->minMatches > numPositions)
			chastise ("--filter can't require more matches (%d) than seed (%d)\n",
					  lzParams->minMatches, numPositions);
		}

	pos_table_set_stat (wordWeight, lzParams->hitSeed->weight);
	pos_table_set_stat (wordSpace,  1L << lzParams->hitSeed->weight);

	lzParams->hitSeed->withTrans = lzParams->withTrans;
	seed_search_set_stat (withTrans,         lzParams->hitSeed->withTrans);
	seed_search_set_stat (minMatches,        lzParams->minMatches);
	seed_search_set_stat (maxTransversions,  lzParams->maxTransversions);
	seed_search_set_stat (filterCaresOnly,   lzParams->filterCaresOnly);

	if (debug >= 90)
		{
		print_seeds (currParams->outputFile, lzParams->hitSeed);
		printf ("%s\n", seed_pattern (lzParams->hitSeed));
		}

	// set span for twins (complicated by the fact that we let the user set it
	// in terms of gap length before we know the seed length)

	if ((lzParams->noHitFiltering) && (twinsYes))
		chastise ("--rawhits can't be used with --twins\n");

	if ((lzParams->noHitFiltering) && (lzParams->gfExtend != gfexNoExtend))
		chastise ("--rawhits can't be used with --gfextend\n");

	if (twinsYes)
		{
		if (minGap <= -lzParams->hitSeed->length)
			chastise ("minGap for twins (%d) must be greater than negative of seed length (%d)\n",
			          minGap, -lzParams->hitSeed->length);
		if (maxGap < minGap)
			chastise ("maxGap for twins (%d) can't be less than min gap (%d)\n",
			          maxGap, minGap);

		lzParams->twinMinSpan = 2 * lzParams->hitSeed->length + minGap;
		lzParams->twinMaxSpan = 2 * lzParams->hitSeed->length + maxGap;
		}
	else
		lzParams->twinMinSpan = 0;  // gaps not used

#ifndef noSeedHitQueue
	if (twinsYes)
		{
		seedHitQueueSize    = lzParams->seedHitQueueSize;
		seedHitQueueColumns = lzParams->twinMaxSpan - lzParams->hitSeed->length;
		}
	else
		{
		seedHitQueueSize    = 0;
		seedHitQueueColumns = -1;
		}
#endif // not noSeedHitQueue
	}

//----------
//
// print_params--
//	Dump some of the user-specified parameters to a file (for debugging).
//
//----------
//
// Arguments:
//	FILE*		f:		The file to print to.
//	control*	params:	Control data to print (some of).
//
// Returns:
//	(nothing)
//
//----------

static void print_params
   (FILE*		f,
	control*	params)
	{
	if (params->seq1 != NULL) fprintf (f, "seq 1: %s\n", params->seq1->fileName);
	if (params->seq2 != NULL) fprintf (f, "seq 2: %s\n", params->seq2->fileName);
	if (params->selfCompare)   fprintf (f, "--self\n");
	print_score_matrix (f, params->scoring, true);
	if      (params->whichStrand > 0) fprintf (f, "--strand=both\n");
	else if (params->whichStrand < 0) fprintf (f, "--strand=minus\n");
	else                               fprintf (f, "--strand=plus\n");
	fprintf (f, "--step=%u\n", params->step);
	fprintf (f, "--seed=%s\n", seed_pattern(params->hitSeed));
	if (params->gfExtend == gfexXDrop) fprintf (f, "--gfextend\n");
	fprintf (f, "--hspthresh=%s\n",    score_thresh_to_string (&params->hspThreshold));
	fprintf (f, "--gappedthresh=%s\n", score_thresh_to_string (&params->gappedThreshold));
	fprintf (f, "--xDrop=" scoreFmtSimple "\n", params->xDrop);
	fprintf (f, "--yDrop=" scoreFmtSimple "\n", params->yDrop);
	fprintf (f, "%s\n", (params->entropicHsp)? "--entropy" : "--noentropy");
	if (params->minMatches >= 0)
		{
		char* qualifier = (params->filterCaresOnly)? "cares:" : "";
		if (params->maxTransversions < 0)
			fprintf (f, "--filter=%s%d\n",    qualifier, params->minMatches);
		else
			fprintf (f, "--filter=%s%d,%d\n", qualifier, params->minMatches, params->maxTransversions);
		}
	if (params->twinMinSpan > 0) fprintf (f, "--twins=%d..%d\n", params->twinMinSpan-2*params->hitSeed->length, params->twinMaxSpan-2*params->hitSeed->length);
	if (params->innerThreshold > 0) fprintf (f, "--innerthresh=" scoreFmtSimple "\n", params->innerThreshold);
	if (params->tracebackMem != defaultParams.tracebackMem) fprintf (f, "--traceback=%d\n", params->tracebackMem);
	}

//----------
//
// read_control_file_by_name, read_control_file--
//	Read control data from a file (see format description below).
//
//----------
//
// Arguments:
//	FILE*		f:		(read_control_file only) The file that control data is
//						.. to be read from.  This should already be open for
//						.. text read.
//	char*		name:	The name of the file that control data is to be read
//						.. from.  For read_control_file this is only used for
//						.. reporting problems to the user (and may be NULL).
//	control*	params:	Control data to fill in.
//
// Returns:
//	(nothing)
//
//----------
//
// Control Data File Format
// ========================
//
// Here's an example:
//
//		min_identity  = 25.0%	# 25th percentile
//		max_identity  = 75.0%	# 75th percentile
//		hsp_threshold = 3000
//
// The control data consists of name-value settings.  Valid names are as
// follows:
//
//		min_identity		The range of sequence identity upon which inference
//		max_identity		.. is based.  Only alignment blocks within this
//							.. range contribute to inference.  If the value ends
//							.. with a percent sign, the range is a percentile of
//							.. the values found in the overall alignment (other-
//							.. wise it is a fixed percentage.
//
//		min_coverage		The range of query coverage upon which inference
//		max_coverage		.. is based.  Only alignment blocks within this
//							.. range (as a percentage of the query sequence)
//							.. contribute to inference.
//
//		min_continuity		The range of query continuity upon which inference
//		max_continuity		.. is based.  Only alignment blocks within this
//							.. range (as a percentage of the query sequence)
//							.. contribute to inference.
//
//		inference_scale		The value for the largest substitution score (i.e.
//							.. the score for the best match).  All other scores
//							.. are scaled by the same factor.  If this is an
//							.. integer (i.e. has no decimal point), then all
//							.. scores will be rounded to an integer as well.
//
//		hsp_threshold		These correspond to the command line --hspthresh
//		gapped_threshold	.. and --gappedthresh options (also known as K and L
//							.. in BLASTZ lingo).  They can be specified as
//							.. mulitples of the scale, e.g.
//							..    hsp_threshold = 20*inference_scale
//							.. Further, the gapped threshold can be specified
//							.. as a multiple of the hsp threshold, e.g.
//							..    gapped_threshold = 1.2*hsp_threshold
//
//		max_sub_iterations	Limits on the number of iterations that will be
//		max_gap_iterations	.. performed.  For example,
//							..    max_sub_iterations = 1
//							..    max_gap_iterations = 0
//							.. will just do one pass and create a substitution
//							.. scoring matrix.
//
//		gap_open_penalty	Correspond to the command line --gap=<[open,]extend>
//		gap_extend_penalty	.. option (also known as O and E in BLASTZ lingo).
//							.. These are the values used for the first iteration
//							.. of gap-scoring inference.  They can be specified
//							.. as mulitples of the scale, and the extend penalty
//							.. can be a multiple of the open penalty.
//
//		step				Corresponds to the command line --step option (also
//							.. known as Z in BLASTZ lingo).
//
//		entropy				Corresponds to the command line --entropy option
//							.. (also known as P in BLASTZ lingo).  Legal values
//							.. are "on" or "off".
//
//----------

static void read_control_file_by_name
   (char*		name,
	control*	params)
	{
	FILE*		f;

	if (name == NULL)
		suicide ("can't open NULL file in read_control_file_by_name()");

	f = fopen_or_die (name, "rt");
	read_control_file (f, name, params);
	fclose_if_valid (f);
	}

static void read_control_file
   (FILE*		f,
	char*		_name,
	control*	params)
	{
	char		line[201];
	char*		name = _name;
	int			lineNum, len, valLen, missingEol;
	char*		waffle;
	char*		valString = NULL;
	int			idIsPercentile = -1;
	int			haveMinIdentity   = false;
	int			haveMaxIdentity   = false;
	int			haveMinCoverage   = false;
	int			haveMaxCoverage   = false;
	int			haveMinContinuity = false;
	int			haveMaxContinuity = false;
	int			tempInt;

	if (name == NULL)
		name = "(unnamed file)";

	//////////
	// read assignments
	//////////

	lineNum    = 0;
	missingEol = false;

	while (fgets (line, sizeof(line), f) != NULL)
		{
		lineNum++;

		// check for lines getting split by fgets (the final line in the file
		// might not have a newline, but no internal lines can be that way)

		if (missingEol)
			goto line_too_long;

		len = strlen (line);
		if (len == 0) continue;
		missingEol = (line[len-1] != '\n');

		// trim blanks, end of line, and comments, and ignore blank lines

		if (line[len-1] == '\n') line[--len] = 0;

		waffle = strchr (line, '#');
		if (waffle != NULL) *waffle = 0;

		trim_string (line);
		if (line[0] == 0) continue;

		// separate the value from the assignment

		valString = strchr (line, '=');
		if (valString == NULL)
			goto invalid_line;

		*(valString++) = 0;
		trim_string (line);
		trim_string (valString);
		valLen = strlen (valString);
		if (valLen == 0)
			goto empty_assignment;

		// parse the assignment

		// inference_scale

		if (strcmp (line, "inference_scale") == 0)
			{
#if (scoreType != 'I')
			int  v;
			char extra;
#endif

			if (strcmp (valString, "none") == 0)
				{
				params->ic.inferScale = 0;
				params->ic.writeAsInt = false;
				}
			else
				{
				params->ic.inferScale = string_to_score (valString);

#if (scoreType == 'I')
				params->ic.writeAsInt = false; // (no need to force it to an int)
#else
				params->ic.writeAsInt = (sscanf (valString, "%d%c", &v, &extra) == 1);
#endif
				}
			continue;
			}

		// hsp_threshold, gapped_threshold
		// $$$ the use of expressions involving other settings should be
		// $$$ .. expanded and generalized

		if (strcmp (line, "hsp_threshold") == 0)
			{
			params->ic.hspThresholdIsRatio = ratioNone;
			if (strcmp_prefix (valString, "top") == 0)
				params->hspThreshold = string_to_score_thresh (valString);
			else if (strcmp_suffix (valString,"*inference_scale") == 0)
				{
				valString[valLen-strlen("*inference_scale")] = 0;
				params->hspThreshold.t = 'S';
				params->hspThreshold.s = string_to_double (valString);
				if (params->ic.inferScale > 0)
					params->hspThreshold.s *= params->ic.inferScale;
				else
					params->ic.hspThresholdIsRatio = ratioMaxSubScore;
				}
			else if (strcmp_suffix (valString,"*worst_substitution") == 0)
				{
				valString[valLen-strlen("*worst_substitution")] = 0;
				params->hspThreshold.t = 'S';
				params->hspThreshold.s = string_to_double (valString);
				params->ic.hspThresholdIsRatio = ratioMinSubScore;
				}
			else
				{
				params->hspThreshold.t = 'S';
				params->hspThreshold.s = string_to_score (valString);
				}
			continue;
			}

		if (strcmp (line, "gapped_threshold") == 0)
			{
			params->ic.gappedThresholdIsRatio = ratioNone;
			if (strcmp_prefix (valString, "top") == 0)
				params->gappedThreshold = string_to_score_thresh (valString);
			else if (strcmp_suffix (valString,"*inference_scale") == 0)
				{
				valString[valLen-strlen("*inference_scale")] = 0;
				params->gappedThreshold.t = 'S';
				params->gappedThreshold.s = string_to_double (valString);
				if (params->ic.inferScale > 0)
					params->gappedThreshold.s *= params->ic.inferScale;
				else
					params->ic.gappedThresholdIsRatio = ratioMaxSubScore;
				}
			else if (strcmp_suffix (valString,"*worst_substitution") == 0)
				{
				valString[valLen-strlen("*worst_substitution")] = 0;
				params->gappedThreshold.t = 'S';
				params->gappedThreshold.s = string_to_double (valString);
				params->ic.gappedThresholdIsRatio = ratioMinSubScore;
				}
			else if (strcmp (valString, "hsp_threshold") == 0)
				params->gappedThreshold = params->hspThreshold;
			else
				{
				params->gappedThreshold.t = 'S';
				params->gappedThreshold.s = string_to_score (valString);
				}
			continue;
			}

		// gap_open_penalty, gap_extend_penalty

		if (strcmp (line, "gap_open_penalty") == 0)
			{
			params->ic.gapOpenIsRatio = ratioNone;
			if (strcmp_suffix (valString,"*inference_scale") == 0)
				{
				valString[valLen-strlen("*inference_scale")] = 0;
				params->scoring->gapOpen = string_to_double (valString);
				if (params->ic.inferScale > 0)
					params->scoring->gapOpen *= params->ic.inferScale;
				else
					params->ic.gapOpenIsRatio = ratioMaxSubScore;
				}
			else if (strcmp_suffix (valString,"*worst_substitution") == 0)
				{
				valString[valLen-strlen("*worst_substitution")] = 0;
				params->scoring->gapOpen = string_to_double (valString);
				params->ic.gapOpenIsRatio = ratioMinSubScore;
				}
			else
				params->scoring->gapOpen = string_to_score (valString);
			continue;
			}

		if (strcmp (line, "gap_extend_penalty") == 0)
			{
			params->ic.gapExtendIsRatio = ratioNone;
			if (strcmp_suffix (valString,"*inference_scale") == 0)
				{
				valString[valLen-strlen("*inference_scale")] = 0;
				params->scoring->gapExtend = string_to_double (valString);
				if (params->ic.inferScale > 0)
					params->scoring->gapExtend *= params->ic.inferScale;
				else
					params->ic.gapExtendIsRatio = ratioMaxSubScore;
				}
			else if (strcmp_suffix (valString,"*worst_substitution") == 0)
				{
				valString[valLen-strlen("*worst_substitution")] = 0;
				params->scoring->gapExtend = string_to_double (valString);
				params->ic.gapExtendIsRatio = ratioMinSubScore;
				}
			else if (strcmp_suffix (valString,"*gap_open_penalty") == 0)
				{
				valString[valLen-strlen("*gap_open_penalty")] = 0;
				params->scoring->gapExtend = string_to_double (valString)
				                           * params->scoring->gapOpen;
				params->ic.gapExtendIsRatio = params->ic.gapOpenIsRatio;
				}
			else
				params->scoring->gapExtend = string_to_score (valString);
			continue;
			}

		// entropy

		if (strcmp (line, "entropy") == 0)
			{
			if (strcmp (valString, "on") == 0)
				params->entropicHsp = true;
			else if (strcmp (valString, "off") == 0)
				params->entropicHsp = false;
			else
				goto on_off_mismatch;
			continue;
			}

		// max_sub_iterations, max_gap_iterations

		if (strcmp (line, "max_sub_iterations") == 0)
			{
			params->ic.subIterations = string_to_int (valString);
			continue;
			}

		if (strcmp (line, "max_gap_iterations") == 0)
			{
			params->ic.gapIterations = string_to_int (valString);
			continue;
			}

		// step

		if (strcmp (line, "step") == 0)
			{
			tempInt = string_to_int (valString);
			if (tempInt <= 0)
				goto bad_step;
			params->step = tempInt;
			continue;
			}

		// min_identity, max_identity
		//
		// Ideally the user will set both of these;  however, if one is set but
		// not the other, we peg the other to the edge of the range;  further,
		// we validate that the user specifies percentile or non-percentile the
		// same for both

		if ((strcmp (line, "min_identity") == 0)
		 || (strcmp (line, "max_identity") == 0))
			{
			if (valString[valLen-1] == '%')
				{
				valString[--valLen] = 0;
				if (idIsPercentile == false)
					goto percentile_mismatch;
				if (idIsPercentile == -1)
					params->ic.idIsPercentile = idIsPercentile = true;
				}
			else if (idIsPercentile == true)
				goto percentile_mismatch;
			else if (idIsPercentile == -1)
				params->ic.idIsPercentile = idIsPercentile = false;

			if (strcmp (line, "min_identity") == 0)
				{
				params->minIdentity = string_to_double (valString) / 100;
				haveMinIdentity = true;
				if (!haveMaxIdentity) params->maxIdentity = 1.0;
				}
			else
				{
				params->maxIdentity = string_to_double (valString) / 100;
				haveMaxIdentity = true;
				if (!haveMinIdentity) params->minIdentity = 0.0;
				}
			continue;
			}

		// min_coverage, max_coverage
		//
		// Ideally the user will set both of these;  however, if one is set but
		// not the other, we peg the other to the edge of the range

		if ((strcmp (line, "min_coverage") == 0)
		 || (strcmp (line, "max_coverage") == 0))
			{
			if (strcmp (line, "min_coverage") == 0)
				{
				params->minCoverage = string_to_double (valString) / 100;
				haveMinCoverage = true;
				if (!haveMaxCoverage) params->maxCoverage = 1.0;
				}
			else
				{
				params->maxCoverage = string_to_double (valString) / 100;
				haveMaxCoverage = true;
				if (!haveMinCoverage) params->minCoverage = 0.0;
				}
			continue;
			}

		// min_continuity, max_continuity
		//
		// Ideally the user will set both of these;  however, if one is set but
		// not the other, we peg the other to the edge of the range

		if ((strcmp (line, "min_continuity") == 0)
		 || (strcmp (line, "max_continuity") == 0))
			{
			if (strcmp (line, "min_continuity") == 0)
				{
				params->minContinuity = string_to_double (valString) / 100;
				haveMinContinuity = true;
				if (!haveMaxContinuity) params->maxContinuity = 1.0;
				}
			else
				{
				params->maxContinuity = string_to_double (valString) / 100;
				haveMaxContinuity = true;
				if (!haveMinContinuity) params->minContinuity = 0.0;
				}
			continue;
			}

		// min_match_count

		if (strcmp (line, "min_match_count") == 0)
			{
			params->minMatchCount = string_to_int (valString);
			continue;
			}

		goto unknown_assignment;
		}

	return;

	//////////
	// failure exits
	//////////

line_too_long:
	suicidef ("line is too long (%s: line %d)", name, lineNum-1);

invalid_line:
	suicidef ("invalid line (%s: line %d) %s", name, lineNum, line);

empty_assignment:
	suicidef ("empty value in assignment (%s: line %d) %s=",
			  name, lineNum, line);

percentile_mismatch:
	suicidef ("assignment of identity percentile mismatches earlier setting"
	          " (%s: line %d) %s=",
	          name, lineNum, line);

on_off_mismatch:
	suicidef ("invalid on/off in assignment (%s: line %d) %s=%s",
	          name, lineNum, line, valString);

unknown_assignment:
	suicidef ("invalid name in assignment (%s: line %d) %s=%s",
	          name, lineNum, line, valString);

bad_step:
	suicidef ("invalid value for step (%s: line %d) %s=%s",
	          name, lineNum, line, valString);
	}

//----------
//
// print_options--
//	Print some of the command line options in the output file.
//
//----------
//
// Arguments:
//	(none)
//
// Returns:
//	(nothing)
//
//----------

static void print_options
   (void)
	{
	print_generic (currParams->outputFile,
	               "seed=%s%s",
	               seed_pattern(currParams->hitSeed),
	               (currParams->hitSeed->withTrans == 0)?"":
	               (currParams->hitSeed->withTrans == 1)?" w/transition"
	                                                    :" w/2 transitions");
	//print_generic (currParams->outputFile, "--hspthresh=" scoreFmtSimple, currParams->hspThreshold);
	//print_generic (currParams->outputFile, "--gappedthresh=" scoreFmtSimple, currParams->gappedThreshold);
	//print_generic (currParams->outputFile, "--xDrop=" scoreFmtSimple, currParams->xDrop);
	//print_generic (currParams->outputFile, "--yDrop=" scoreFmtSimple, currParams->yDrop);
	//print_generic (currParams->outputFile, "%s", (currParams->entropicHsp)? "--entropy" : "--noentropy");
	//if (currParams->minMatches >= 0) print_generic (currParams->outputFile, 'z', "--filter=%d,%d", currParams->minMatches, currParams->maxTransversions);
	//if (currParams->twinMinSpan > 0) print_generic (currParams->outputFile, "twins=%d..%d", currParams->twinMinSpan-2*currParams->hitSeed->length, currParams->twinMaxSpan-2*currParams->hitSeed->length);
	//                        else print_generic (currParams->outputFile, "notwins");
	print_generic (currParams->outputFile, "step=%u", currParams->step);
	}

//----------
//
// name_spec_is_quantum--
//	Determine if a sequence name specifier is describing a quantum sequence.
//
// It is deemed to be a quantum sequence if either the filename ends in ".qdna"
// or the action list contains the "quantum" action.
//	Dump the nucleotides (from each sequence) for a gap-free alignment.
//
//----------
//
// Arguments:
//	char*	spec:	The sequence name specifier.  This is of the form
//					.. <filename>[<action_lists>], where each action list is a
//					.. comma-separated list inside square brackets.
//
// Returns:
//	(nothing)
//
//----------

static int name_spec_is_quantum
   (char*	spec)
	{
	char*	nameEnd, *where;
	char	before, after;

	if (spec == NULL) return false;

	// see if file name ends with ".qdna"

	nameEnd = strchr(spec,'[');
	if (nameEnd == NULL)
		return  (strcmp_suffix(spec,".qdna") == 0);

	if (strncmp_suffix(spec,".qdna",nameEnd-spec) == 0) return true;

	// see if action lists contain the keyword "quantum"

	where = strstr(nameEnd,"quantum");		// (nota bene: nameEnd != NULL)
	if (where != NULL)
		{
		before = where[-1];					// (nota bene: where > nameEnd)
		after  = where[strlen("quantum")];
		if (((before == '[') || (before == ','))
		 && ((after  == ']') || (after  == ',') || (after == '=')))
			return true;
		}

	return false;
	}

//----------
//
// lastz_zero_stats--
//	Clear the statistics for this module.
//
//----------
//
// Arguments:
//	(none)
//
// Returns:
//	(nothing)
//
//----------

void lastz_zero_stats
   (void)
	{
#ifdef collect_stats

	// set 'em en masse to zero

	memset (&lastzStats, 0, sizeof(lastzStats));

	// set any values that might be floating point to zero (fp bit pattern for
	// zero may not be all-bits-zero)

	// (none to set, yet)

#endif // collect_stats
	}

//----------
//
// lastz_show_stats_before, lastz_show_stats--
//	Show the statistics that have been collected for this module.
//
//----------
//
// Arguments:
//	FILE*		f:	The file to print the stats to.
//
// Returns:
//	(nothing)
//
//----------

static void lastz_show_stats_before
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats
	if (f == NULL) return;
	fprintf (f, "-------------------\n");
	fprintf (f, "     target length: %s\n",     commatize (lzParams.seq1->len));
	if (lzParams.seq2 != NULL)
		fprintf (f, "      query length: %s\n", commatize (lzParams.seq2->len));
	fprintf (f, "         step size: %u\n",     currParams->step);
	fprintf (f, "-------------------\n");
#endif // collect_stats
	}


static void lastz_show_stats
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats
	if (f == NULL) return;
	fprintf (f, "          run time: %.3f seconds\n", lastzStats.runTime);
	fprintf (f, "-------------------\n");
#endif // collect_stats
	}

