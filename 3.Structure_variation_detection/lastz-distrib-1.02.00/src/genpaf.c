//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: genpaf.c
//
//----------
//
// genpaf--
//	Support for printing alignments in "GENeral Pairwise Alignment Format".
//
// genpaf format is non-standard.  It prints each alignment block on a single
// line, and the calling program can specify which fields are printed, and in
// what order.  It is best suited for situations in which the alignment file
// will be processed by some other program.
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
#include <string.h>				// standard C string stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff
#include "diag_hash.h"			// diagonals hashing stuff
#include "identity_dist.h"		// identity distribution stuff
#include "coverage_dist.h"		// query coverage distribution stuff
#include "continuity_dist.h"	// query continuity distribution stuff
#include "cigar.h"				// cigar alignment format stuff

#define  genpaf_owner			// (make this the owner of its globals)
#include "genpaf.h"				// interface to this module

//----------
//
// prototypes for private functions
//
//----------

static char* extract_key_info (char* keys, char desiredKey);

//----------
//
// print_genpaf_job_header--
//	Print genpaf format job header.
//
//----------
//
// Arguments:
//	FILE*	f:		The file to print to.
//	char*	keys:	A list of the fields we'll be printing (genpafXXX values).
//
// Returns:
//	(nothing)
//
//----------


void print_genpaf_job_header
   (FILE*	f,
	char*	keys)
	{
	char*	k;
	int		tabCh;

	if (keys == NULL) return;

	// print headers for the desired fields
	
	tabCh = '#';
	for (k=keys ; (*k!=0)&&(*k!=genpafInfoSeparator); k++)
		{
		if      (tabCh == '#')       { fprintf (f, "#");  tabCh = '\t'; }
		else if (tabCh == 0)         { fprintf (f, "#");  tabCh = '\t'; }
		else if (*k == genpafCR)     {                    tabCh = '\t'; }
		else if (*k == genpafMarker) {                    tabCh = '\t'; }
		                        else { fprintf (f, "\t");               }
		switch (*k)
			{
			case genpafCR:               fprintf (f, "\n");  tabCh = '#';   break;
			case genpafMarker:           fprintf (f, "~");   tabCh = 0;     break;
			case genpafNA:                                                  break;
			case genpafName1:            fprintf (f, "name1");              break;
			case genpafStrand1:          fprintf (f, "strand1");            break;
			case genpafSize1:            fprintf (f, "size1");              break;
			case genpafStart1:           fprintf (f, "start1");             break;
			case genpafStart1Zero:       fprintf (f, "zstart1");            break;
			case genpafStart1DotPlot:    fprintf (f, "start1");             break;
			case genpafEnd1:             fprintf (f, "end1");               break;
			case genpafEnd1DotPlot:      fprintf (f, "end1");               break;
			case genpafLength1:          fprintf (f, "length1");            break;
			case genpafAlign1:           fprintf (f, "align1");             break;
			case genpafText1:            fprintf (f, "text1");              break;
			case genpafName2:            fprintf (f, "name2");              break;
			case genpafStrand2:          fprintf (f, "strand2");            break;
			case genpafSize2:            fprintf (f, "size2");              break;
			case genpafStart2:           fprintf (f, "start2");             break;
			case genpafStart2Zero:       fprintf (f, "zstart2");            break;
			case genpafStart2OnPlus:     fprintf (f, "start2+");            break;
			case genpafStart2ZeroOnPlus: fprintf (f, "zstart2+");           break;
			case genpafStart2DotPlot:    fprintf (f, "start2");             break;
			case genpafEnd2:             fprintf (f, "end2");               break;
			case genpafEnd2OnPlus:       fprintf (f, "end2+");              break;
			case genpafEnd2DotPlot:      fprintf (f, "end2");               break;
			case genpafLength2:          fprintf (f, "length2");            break;
			case genpafAlign2:           fprintf (f, "align2");             break;
			case genpafText2:            fprintf (f, "text2");              break;
			case genpafMatch:            fprintf (f, "nmatch");             break;
			case genpafMismatch:         fprintf (f, "nmismatch");          break;
			case genpafSeparateGaps:     fprintf (f, "ngap");               break;
			case genpafGapColumns:       fprintf (f, "cgap");               break;
			case genpafTextDiff:         fprintf (f, "diff");               break;
			case genpafCigar:            fprintf (f, "cigar");              break;
			case genpafCigarLower:       fprintf (f, "cigar-");             break;
			case genpafCigarX:           fprintf (f, "cigarx");             break;
			case genpafCigarXLower:      fprintf (f, "cigarx-");            break;
			case genpafDiagonal:         fprintf (f, "diagonal");           break;
			case genpafShingle:          fprintf (f, "shingle");            break;
			case genpafScore:            fprintf (f, "score");              break;
			case genpafIdentity:         fprintf (f, "identity\tidPct");    break;
			case genpafCoverage:         fprintf (f, "coverage\tcovPct");   break;
			case genpafContinuity:       fprintf (f, "continuity\tconPct"); break;
			case genpafGapRate:          fprintf (f, "gaprate\tgapPct");    break;
			default:                                                        break;
			}
		}
	fprintf (f, "\n");
	}

//----------
//
// print_genpaf_job_footer--
//	Print genpaf format job footer.
//
//----------

void print_genpaf_job_footer
   (arg_dont_complain(FILE* f))
	{
	// (do nothing)
	}

//----------
//
// print_genpaf_header--
//	Print genpaf format query header.
//
//----------

void print_genpaf_header
   (arg_dont_complain(FILE* f),
	arg_dont_complain(seq*  seq1),
	arg_dont_complain(seq*  seq2))
	{
	// (do nothing)
	}

//----------
//
// print_genpaf_align_list--
//	Print a list of gapped alignments in genpaf format.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	alignel*	alignList:	The list of alignments to print.
//	seq*		seq1:		One sequence.
//	seq*		seq2:		Another sequence.
//	char*		keys:		A list of the fields to print (genpafXXX values).
//
// Returns:
//	(nothing)
//
//----------

void print_genpaf_align_list
   (FILE*		f,
	alignel*	alignList,
	seq*		seq1,
	seq*		seq2,
	char*		keys)
	{
	alignel*	a;
	int			computeId, computeCov, computeCon, computeGap;
	unspos		idNumer, idDenom;
	unspos		covNumer, covDenom;
	unspos		conNumer, conDenom;
	unspos		gapNumer, gapDenom;

	idNumer = idDenom = covNumer = covDenom = 0;

	computeId  = (strchr (keys, genpafIdentity) != NULL);
	computeCov = (strchr (keys, genpafCoverage) != NULL);
	computeCon = (strchr (keys, genpafContinuity)  != NULL);
	computeGap = (strchr (keys, genpafGapRate)  != NULL);

	for (a=alignList ; a!=NULL ; a=a->next)
		{
		if (computeId)  alignment_identity   (seq1, seq2, a, &idNumer,  &idDenom);
		if (computeCov) alignment_coverage   (seq1, seq2, a, &covNumer, &covDenom);
		if (computeCon) alignment_continuity (            a, &conNumer, &conDenom);
		if (computeGap) alignment_gap_rate   (            a, &gapNumer, &gapDenom);

		print_genpaf_align (f,
		                    seq1, a->beg1-1, a->end1,
		                    seq2, a->beg2-1, a->end2,
		                    a->script, a->s,
		                    keys,
		                    idNumer,  idDenom,
		                    covNumer, covDenom,
		                    conNumer, conDenom,
		                    gapNumer, gapDenom);
		}
	}

//----------
//
// print_genpaf_align_list_segments--
//	Print a list of gapped alignments in genpaf format, splitting them into
//	ungapped segments.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	alignel*	alignList:	The list of alignments to print.
//	seq*		seq1:		One sequence.
//	seq*		seq2:		Another sequence.
//	char*		keys:		A list of the fields to print (genpafXXX values).
//
// Returns:
//	(nothing)
//
//----------

void print_genpaf_align_list_segments
   (FILE*		f,
	alignel*	alignList,
	seq*		seq1,
	seq*		seq2,
	char*		keys)
	{
	alignel*	a;
	unspos		beg1, end1, beg2, end2;
	unspos		height, width, i, j, prevI, prevJ, run;
	u32			opIx;

	for (a=alignList ; a!=NULL ; a=a->next)
		{
		beg1   = a->beg1;
		end1   = a->end1;
		beg2   = a->beg2;
		end2   = a->end2;
		height = end1 - beg1 + 1;
		width  = end2 - beg2 + 1;

		// print the alignment's segments

		opIx = 0;
		for (i=j=0 ; (i< height)||(j<width) ; )
			{
			prevI = i;  prevJ = j;
			run = edit_script_run_of_subs (a->script, &opIx);
			i += run; j += run;
			if ((i < height) || (j < width))
				edit_script_indel_len (a->script, &opIx, &i, &j);

			print_genpaf_match (f,
			                    seq1, beg1-1+prevI,
			                    seq2, beg2-1+prevJ,
			                    run, /*score*/ 0,
			                    keys);
			}
		}

	}

//----------
//
// print_genpaf_align--
//	Print a single gapped alignment in genpaf format.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	seq*		seq1:		One sequence.
//	unspos		beg1, end1:	Range of positions in sequence 1 (origin 0).
//	seq*		seq2:		Another sequence.
//	unspos		beg2, end2:	Range of positions in sequence 2 (origin 0).
//	editscript*	script:		The script describing the path the alignment
//							.. takes in the DP matrix.
//	score		s:			The alignment's score.
//	char*		keys:		A list of the fields to print (genpafXXX values).
//	unspos		idNumer:	Identity
//	unspos		idDenom:
//	unspos		covNumer:	Coverage
//	unspos		covDenom:
//	unspos		gapNumer:	Gap Rate
//	unspos		gapDenom:
//
// Returns:
//	(nothing)
//
//----------

static char* rcfSuffix[4] = { "", "~", "~", "" };


static char diff_char (u8 p, u8 q, char* textDiffInfo);
static char diff_char (u8 p, u8 q, char* textDiffInfo)
	{
	s8		b1,  b2;
	char	c;

	b1 = nuc_to_bits[p];
	b2 = nuc_to_bits[q];

	if      ((b1 < 0) || (b2 < 0))                               c = textDiffInfo[genpafTDInfoOther];
	else if (b1 == b2)                                           c = textDiffInfo[genpafTDInfoMatch];
	else if (bits_to_pur_pyr[(u8)b1] == bits_to_pur_pyr[(u8)b2]) c = textDiffInfo[genpafTDInfoTransition];
	else                                                         c = textDiffInfo[genpafTDInfoTransversion];

	return c;
	}


void print_genpaf_align
   (FILE*			f,
	seq*			seq1,
	unspos			beg1,
	unspos			end1,
	seq*			seq2,
	unspos			beg2,
	unspos			end2,
	editscript*		script,
	score			s,
	char*			keys,
	unspos			idNumer,
	unspos			idDenom,
	unspos			covNumer,
	unspos			covDenom,
	unspos			conNumer,
	unspos			conDenom,
	unspos			gapNumer,
	unspos			gapDenom)
	{
	char*			k;
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	unspos			height, width, i, j, prevI, prevJ, run;
	u32				opIx;
	u8*				p, *q;
	unspos			ix, len1, len2;
	char*			name1, *name2, *suff1, *suff2;
	unspos			offset1, offset2, start1, start2;
	unspos			dotStart1, dotStart2, dotEnd1, dotEnd2;
	unspos			seq1Len, seq2Len, seq1True, seq2True;
	char			strand1, strand2;
	unspos			startI, startJ;
	int				tabCh;
	char*			textDiffInfo;
	sgnpos			diag, diagSE, diagNW;
	unspos			numGaps;

	beg1++; // (internally, we want origin 1, inclusive)
	beg2++;

	len1 = height = end1 - beg1 + 1;
	len2 = width  = end2 - beg2 + 1;

	//////////
	// figure out position offsets and names
	//////////

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		name1 = (seq1->useFullNames)? seq1->header : seq1->shortHeader;
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";
		offset1  = 0;
		seq1Len  = seq1->len;
		seq1True = seq1->trueLen;
		}
	else					// sequence 1 is partitioned
	 	{
		part = lookup_partition (seq1, beg1-1);
		name1    = &sp1->pool[part->header];
		offset1  = part->sepPos + 1;
		seq1Len  = (part+1)->sepPos - offset1;
		seq1True = part->trueLen;
		}

	if (sp2->p == NULL)		// sequence 2 is not partitioned
		{
		name2 = (seq2->useFullNames)? seq2->header : seq2->shortHeader;
		if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";
		offset2  = 0;
		seq2Len  = seq2->len;
		seq2True = seq2->trueLen;
		}
	else					// sequence 2 is partitioned
	 	{
		part = lookup_partition (seq2, beg2-1);
		name2    = &sp2->pool[part->header];
		offset2  = part->sepPos + 1;
		seq2Len  = (part+1)->sepPos - offset2;
		seq2True = part->trueLen;
		}

	//////////
	// figure out strandedness
	//////////

	suff1 = rcfSuffix[seq1->revCompFlags];
	suff2 = rcfSuffix[seq2->revCompFlags];

	if ((seq1->revCompFlags & rcf_rev) == 0)
		{
		start1    = beg1-1 - offset1 + seq1->start;
		dotStart1 = start1 + 1;
		dotEnd1   = dotStart1 + height - 1;
		strand1   = '+';
		}
	else
		{
		start1  = beg1-1 - offset1 + seq1True+2 - (seq1->start + seq1Len);
		dotStart1 = (seq1->start + seq1Len + offset1 - beg1) - 1;
		dotEnd1   = (dotStart1 - height) + 1;
		strand1 = '-';
		}
	if ((seq2->revCompFlags & rcf_rev) == 0)
		{
		start2    = beg2-1 - offset2 + seq2->start;
		dotStart2 = start2 + 1;
		dotEnd2   = dotStart2 + width - 1;
		strand2   = '+';
		}
	else
		{
		start2    = beg2-1 - offset2 + seq2True+2 - (seq2->start + seq2Len);
		dotStart2 = (seq2->start + seq2Len + offset2 - beg2) - 1;
		dotEnd2   = (dotStart2 - width) + 1;
		strand2   = '-';
		}

	// print the desired fields

	textDiffInfo = extract_key_info (keys, genpafTextDiff);
	if (textDiffInfo == NULL) textDiffInfo = genpafTDInfoDefault;

	tabCh = '#';
	for (k=keys ; (*k!=0)&&(*k!=genpafInfoSeparator); k++)
		{
		if ((tabCh == '#')
		 || (tabCh == 0)
		 || (*k == genpafCR)
		 || (*k == genpafMarker))
			tabCh = '\t';
		else
			fprintf (f, "\t");

		switch (*k)
			{
			case genpafCR:
				fprintf (f, "\n");
				tabCh = '#';
				break;
			case genpafMarker:
				fprintf (f, "~");
				tabCh = 0;
				break;
			case genpafNA:
				fprintf (f, "NA");
				break;
			case genpafName1:
				fprintf (f, "%s%s", name1, suff1);
				break;
			case genpafStrand1:
				fprintf (f, "%c", strand1);
				break;
			case genpafSize1:
				fprintf (f, unsposFmt, seq1True);
				break;
			case genpafStart1:
				fprintf (f, unsposFmt, start1);
				break;
			case genpafStart1Zero:
				fprintf (f, unsposFmt, start1-1);
				break;
			case genpafStart1DotPlot:
				fprintf (f, unsposFmt, dotStart1);
				break;
			case genpafEnd1:
				fprintf (f, unsposFmt, start1+len1-1);
				break;
			case genpafEnd1DotPlot:
				fprintf (f, unsposFmt, dotEnd1);
				break;
			case genpafLength1:
				fprintf (f, unsposFmt, height);
				break;
			case genpafAlign1:
			case genpafText1:
				// print aligning path in sequence 1 (non-printables are printed as '*'
				// but such should never be seen unless there is a problem elsewhere)

				opIx = 0;
				for (i=j=0 ; (i<height)||(j<width) ; )
					{
					// handle the next run

					run = edit_script_run_of_subs (script, &opIx);

					p = seq1->v+beg1+i-1;
					q = seq2->v+beg2+j-1;
					for (ix=0 ; ix<run ; ix++)
						{ fprintf (f, "%c", dna_toprint(*p));  p++;  q++; }

					i += run; j += run;

					// handle the next indel

					if ((i < height) || (j < width))
						{
						startI = i;  p = seq1->v+beg1+i-1;
						startJ = j;  q = seq2->v+beg2+j-1;

						edit_script_indel_len (script, &opIx, &i, &j);

						if (i != startI)
							{
							for ( ; startI<i ; startI++)
								{ fprintf (f, "%c", dna_toprint(*p));  p++; }
							}

						if (j != startJ)
							{
							for ( ; startJ<j ; startJ++)
								{ fprintf (f, "-");  q++; }
							}
						}
					}
				break;
			case genpafName2:
				fprintf (f, "%s%s", name2, suff2);
				break;
			case genpafStrand2:
				fprintf (f, "%c", strand2);
				break;
			case genpafSize2:
				fprintf (f, unsposFmt, seq2True);
				break;
			case genpafStart2OnPlus:
				if (strand2 == '-')
					{
					fprintf (f, unsposFmt, seq2True + 2 - start2 - len2);
					break;
					}
				// (fall thru)
			case genpafStart2:
				fprintf (f, unsposFmt, start2);
				break;
			case genpafStart2ZeroOnPlus:
				if (strand2 == '-')
					{
					fprintf (f, unsposFmt, seq2True + 1 - start2 - len2);
					break;
					}
				// (fall thru)
			case genpafStart2Zero:
				fprintf (f, unsposFmt, start2-1);
				break;
			case genpafStart2DotPlot:
				fprintf (f, unsposFmt, dotStart2);
				break;
			case genpafEnd2OnPlus:
				if (strand2 == '-')
					{
					fprintf (f, unsposFmt, seq2True + 1 - start2);
					break;
					}
				// (fall thru)
			case genpafEnd2:
				fprintf (f, unsposFmt, start2+len2-1);
				break;
			case genpafEnd2DotPlot:
				fprintf (f, unsposFmt, dotEnd2);
				break;
			case genpafLength2:
				fprintf (f, unsposFmt, width);
				break;
			case genpafAlign2:
			case genpafText2:
				// print aligning path in sequence 2

				opIx = 0;
				for (i=j=0 ; (i<height)||(j<width) ; )
					{
					// handle the next run

					run = edit_script_run_of_subs (script, &opIx);

					p = seq1->v+beg1+i-1;
					q = seq2->v+beg2+j-1;
					for (ix=0 ; ix<run ; ix++)
						{ fprintf (f, "%c", dna_toprint(*q));  p++;  q++; }

					i += run; j += run;

					// handle the next indel

					if ((i < height) || (j < width))
						{
						startI = i;  p = seq1->v+beg1+i-1;
						startJ = j;  q = seq2->v+beg2+j-1;

						edit_script_indel_len (script, &opIx, &i, &j);

						if (i != startI)
							{
							for ( ; startI<i ; startI++)
								{ fprintf (f, "-");  p++; }
							}

						if (j != startJ)
							{
							for ( ; startJ<j ; startJ++)
								{ fprintf (f, "%c", dna_toprint(*q));  q++; }
							}
						}
					}
				break;
			case genpafMatch:
				fprintf (f, unsposFmt, idNumer);
				break;
			case genpafMismatch:
				fprintf (f, unsposFmt, idDenom - idNumer);
				break;
			case genpafSeparateGaps:
				// count gaps in aligning path;  note that if we have a
				// single path operator that is a gap in both sequences, we
				// only count this as 1 gap

				numGaps = 0;

				opIx = 0;
				for (i=j=0 ; (i<height)||(j<width) ; )
					{
					// handle the next run

					run = edit_script_run_of_subs (script, &opIx);
					i += run; j += run;

					// handle the next indel

					if ((i < height) || (j < width))
						{
						edit_script_indel_len (script, &opIx, &i, &j);
						numGaps++;
						}
					}
				fprintf (f, unsposFmt, numGaps);
				break;
			case genpafGapColumns:
				fprintf (f, unsposFmt, conDenom - conNumer);
				break;
			case genpafTextDiff:
				// print differences between aligning paths

				opIx = 0;
				for (i=j=0 ; (i<height)||(j<width) ; )
					{
					// handle the next run

					run = edit_script_run_of_subs (script, &opIx);

					p = seq1->v+beg1+i-1;
					q = seq2->v+beg2+j-1;
					for (ix=0 ; ix<run ; ix++)
						{ fprintf (f, "%c", diff_char (*p, *q, textDiffInfo));  p++;  q++; }

					i += run; j += run;

					// handle the next indel

					if ((i < height) || (j < width))
						{
						startI = i;  p = seq1->v+beg1+i-1;
						startJ = j;  q = seq2->v+beg2+j-1;

						edit_script_indel_len (script, &opIx, &i, &j);

						if (i != startI)
							{
							for ( ; startI<i ; startI++)
								{
								fprintf (f, "%c", textDiffInfo[genpafTDInfoInsert1]); 
								p++;
								}
							}

						if (j != startJ)
							{
							for ( ; startJ<j ; startJ++)
								{
								fprintf (f, "%c", textDiffInfo[genpafTDInfoInsert2]); 
								q++;
								}
							}
						}
					}
				break;
			case genpafCigar:
			case genpafCigarLower:
				opIx = 0;
				for (i=j=0 ; (i< height)||(j<width) ; )
					{
					run = edit_script_run_of_subs (script, &opIx);
					if (*k == genpafCigar)
						fprintf (f, unsposFmt "M", run);
					else // if (*k == genpafCigarLower)
						fprintf (f, unsposFmt "m", run);
					i += run; j += run;
			
					if ((i < height) || (j < width))
						{
						prevI = i;  prevJ = j;
						edit_script_indel_len (script, &opIx, &i, &j);
						if (i > prevI)
							{
							if (*k == genpafCigar)
								fprintf (f, unsposFmt "D", i - prevI);
							else
								fprintf (f, unsposFmt "d", i - prevI);
							}
						if (j > prevJ)
							{
							if (*k == genpafCigar)
								fprintf (f, unsposFmt "I", j - prevJ);
							else
								fprintf (f, unsposFmt "i", j - prevJ);
							}
						}
					}
				break;
			case genpafCigarX:
			case genpafCigarXLower:
				print_cigar_align (f, seq1, beg1-1, end1, seq2, beg2-1, end2,
				                   script, s,
				                   /* withInfo       */ false,
				                   /* markMismatches */ true,
				                   /* letterAfter    */ true,
				                   /* hideSingles    */ true,
				                   /* lowerCase      */ (*k == genpafCigarXLower),
				                   /* withNewLine    */ false);
				break;
			case genpafDiagonal:
				fprintf (f, sgnposFmt, diagNumber(start1,start2));
				break;
			case genpafShingle:
				diag = diagNumber (start1, start2);
				diagSE = seq1Len - diag;
				diagNW = seq2Len + diag;
				if (diag < 0)
					{
					if ((diagNW < 0)
					 || ((u32) diagNW < seq1Len)) diag = -diagNW;
					                         else diag = 0;
					}
				else if (diag > 0)
					{
					if ((diagSE < 0)
					 || ((u32) diagSE < seq2Len)) diag = diagSE;
					                         else diag = 0;
					}
				if (diag == 0) fprintf (f, "NA");
				          else fprintf (f, sgnposFmt, diag);
				break;
			case genpafScore:
				fprintf (f, scoreFmt, s);
				break;
			case genpafIdentity:
				fprintf (f, unsposSlashFmt, idNumer, idDenom);
				if (idDenom != 0) fprintf (f, "\t%.1f%%", (100.0*idNumer) / idDenom);
				              else fprintf (f, "\tNA");
				break;
			case genpafCoverage:
				fprintf (f, unsposSlashFmt, covNumer, covDenom);
				if (covDenom != 0) fprintf (f, "\t%.1f%%", (100.0*covNumer) / covDenom);
				              else fprintf (f, "\tNA");
				break;
			case genpafContinuity:
				fprintf (f, unsposSlashFmt, conNumer, conDenom);
				if (conDenom != 0) fprintf (f, "\t%.1f%%", (100.0*conNumer) / conDenom);
				              else fprintf (f, "\tNA");
				break;
			case genpafGapRate:
				fprintf (f, unsposSlashFmt, gapNumer, gapDenom);
				if (gapDenom != 0) fprintf (f, "\t%.1f%%", (100.0*gapNumer) / gapDenom);
				              else fprintf (f, "\tNA");
				break;
			default:
				break;
			}
		}
	fprintf (f, "\n");
	}

//----------
//
// print_genpaf_match--
//	Print an hsp in genpaf format.
//
//----------
//
// Arguments:
//	FILE*	f:		The file to print to.
//	seq*	seq1:	One sequence.
//	unspos	pos1:	The position, in seq1, of first character in the match
//					.. (origin-0).
//	seq*	seq2:	Another sequence.
//	unspos	pos1:	The position, in seq2, of first character in the match
//					.. (origin-0).
//	unspos	length:	The number of nucleotides in the HSP.
//	score	s:		The HSP's score.
//	char*	keys:	A list of the fields to print (genpafXXX values).
//
// Returns:
//	(nothing)
//
//----------

void print_genpaf_match
   (FILE*			f,
	seq*			seq1,
	unspos			pos1,
	seq*			seq2,
	unspos			pos2,
	unspos			length,
	score			s,
	char*			keys)
	{
	char*			k;
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	u8*				s1 = seq1->v + pos1;
	u8*				s2 = seq2->v + pos2;
	char*			name1, *name2, *suff1, *suff2;
	unspos			offset1, offset2, start1, start2;
	unspos			dotStart1, dotStart2, dotEnd1, dotEnd2;
	unspos			seq1Len, seq2Len, seq1True, seq2True;
	char			strand1, strand2;
	unspos			ix;
	segment			seg;
	unspos			numer, denom;
	int				tabCh;
	char*			textDiffInfo;
	sgnpos			diag, diagSE, diagNW;

	if (seq1->revCompFlags != rcf_forward)
		suicide ("attempt to print - strand or complement for sequence 1 in print_genpaf_match");

	//////////
	// figure out position offsets and names
	//////////

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		name1 = (seq1->useFullNames)? seq1->header : seq1->shortHeader;
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";
		offset1  = 0;
		seq1Len  = seq1->len;
		seq1True = seq1->trueLen;
		}
	else					// sequence 1 is partitioned
	 	{
		part = lookup_partition (seq1, pos1);
		name1    = &sp1->pool[part->header];
		offset1  = part->sepPos + 1;
		seq1Len  = (part+1)->sepPos - offset1;
		seq1True = part->trueLen;
		}

	if (sp2->p == NULL)		// sequence 2 is not partitioned
		{
		name2 = (seq2->useFullNames)? seq2->header : seq2->shortHeader;
		if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";
		offset2  = 0;
		seq2Len  = seq2->len;
		seq2True = seq2->trueLen;
		}
	else					// sequence 2 is partitioned
	 	{
		part = lookup_partition (seq2, pos2);
		name2    = &sp2->pool[part->header];
		offset2  = part->sepPos + 1;
		seq2Len  = (part+1)->sepPos - offset2;
		seq2True = part->trueLen;
		}

	//////////
	// figure out strandedness
	//////////

	suff1 = rcfSuffix[seq1->revCompFlags];
	suff2 = rcfSuffix[seq2->revCompFlags];

	if ((seq1->revCompFlags & rcf_rev) == 0)
		{
		start1    = pos1 - offset1 + seq1->start;
		dotStart1 = start1;
		dotEnd1   = dotStart1 + length - 1;
		strand1   = '+';
		}
	else
		{
		start1    = pos1 - offset1 + seq1True+2 - (seq1->start + seq1Len);
		dotStart1 = (seq1->start + seq1Len + offset1 - pos1) - 1;
		dotEnd1   = (dotStart1 - length) + 1;
		strand1   = '-';
		}
	if ((seq2->revCompFlags & rcf_rev) == 0)
		{
		start2    = pos2 - offset2 + seq2->start;
		dotStart2 = start2;
		dotEnd2   = (dotStart2 + length) - 1;
		strand2   = '+';
		}
	else
		{
		start2    = pos2 - offset2 + seq2True+2 - (seq2->start + seq2Len);
		dotStart2 = (seq2->start + seq2Len + offset2 - pos2) - 1;
		dotEnd2   = (dotStart2 - length) + 1;
		strand2   = '-';
		}

	// print the desired fields

	textDiffInfo = extract_key_info (keys,genpafTextDiff);
	if (textDiffInfo == NULL) textDiffInfo = genpafTDInfoDefault;

	tabCh = '#';
	for (k=keys ; (*k!=0)&&(*k!=genpafInfoSeparator); k++)
		{
		if ((tabCh == '#')
		 || (tabCh == 0)
		 || (*k == genpafCR)
		 || (*k == genpafMarker))
			tabCh = '\t';
		else
			fprintf (f, "\t");

		switch (*k)
			{
			case genpafCR:
				fprintf (f, "\n");
				tabCh = '#';
				break;
			case genpafMarker:
				fprintf (f, "~");
				tabCh = 0;
				break;
			case genpafNA:
				fprintf (f, "NA");
				break;
			case genpafName1:
				fprintf (f, "%s%s", name1, suff1);
				break;
			case genpafStrand1:
				fprintf (f, "%c", strand1);
				break;
			case genpafSize1:
				fprintf (f, unsposFmt, seq1True);
				break;
			case genpafStart1:
				fprintf (f, unsposFmt, start1);
				break;
			case genpafStart1Zero:
				fprintf (f, unsposFmt, start1-1);
				break;
			case genpafStart1DotPlot:
				fprintf (f, unsposFmt, dotStart1);
				break;
			case genpafEnd1:
				fprintf (f, unsposFmt, start1-1 + length);
				break;
			case genpafEnd1DotPlot:
				fprintf (f, unsposFmt, dotEnd1);
				break;
			case genpafLength1:
				fprintf (f, unsposFmt, length);
				break;
			case genpafAlign1:
			case genpafText1:
				for (ix=0 ; ix<length ; ix++)
					fprintf (f, "%c", dna_toprint(s1[ix]));
				break;
			case genpafName2:
				fprintf (f, "%s%s", name2, suff2);
				break;
			case genpafStrand2:
				fprintf (f, "%c", strand2);
				break;
			case genpafSize2:
				fprintf (f, unsposFmt, seq2True);
				break;
			case genpafStart2OnPlus:
				if (strand2 == '-')
					{
					fprintf (f, unsposFmt, seq2True + 2 - (start2+length));
					break;
					}
				// (fall thru)
			case genpafStart2:
				fprintf (f, unsposFmt, start2);
				break;
			case genpafStart2ZeroOnPlus:
				if (strand2 == '-')
					{
					fprintf (f, unsposFmt, seq2True + 1 - (start2+length));
					break;
					}
				// (fall thru)
			case genpafStart2Zero:
				fprintf (f, unsposFmt, start2-1);
				break;
			case genpafStart2DotPlot:
				fprintf (f, unsposFmt, dotStart2);
				break;
			case genpafEnd2OnPlus:
				if (strand2 == '-')
					{
					fprintf (f, unsposFmt, seq2True + 1 - start2);
					break;
					}
				// (fall thru)
			case genpafEnd2:
				fprintf (f, unsposFmt, start2-1 + length);
				break;
			case genpafEnd2DotPlot:
				fprintf (f, unsposFmt, dotEnd2);
				break;
			case genpafLength2:
				fprintf (f, unsposFmt, length);
				break;
			case genpafAlign2:
			case genpafText2:
				for (ix=0 ; ix<length ; ix++)
					fprintf (f, "%c", dna_toprint(s2[ix]));
				break;
			case genpafMatch:
				segment_identity (seq1, pos1, seq2, pos2, length, &numer, &denom);
				fprintf (f, unsposFmt, numer);
				break;
			case genpafMismatch:
				segment_identity (seq1, pos1, seq2, pos2, length, &numer, &denom);
				fprintf (f, unsposFmt, denom - numer);
				break;
			case genpafSeparateGaps:
			case genpafGapColumns:
				fprintf (f, "0");
				break;
			case genpafTextDiff:
				for (ix=0 ; ix<length ; ix++)
					fprintf (f, "%c", diff_char (s1[ix], s2[ix], textDiffInfo));
				break;
			case genpafCigar:
				fprintf (f, unsposFmt "M", length);
				break;
			case genpafCigarLower:
				fprintf (f, unsposFmt "m", length);
				break;
			case genpafCigarX:
			case genpafCigarXLower:
				print_cigar_match (f, seq1, pos1, seq2, pos2, length,
				                   s,
				                   /* withInfo       */ false,
				                   /* markMismatches */ true,
				                   /* letterAfter    */ true,
				                   /* hideSingles    */ true,
				                   /* lowerCase      */ (*k == genpafCigarXLower),
				                   /* withNewLine    */ false);
				break;
			case genpafDiagonal:
				fprintf (f, sgnposFmt, diagNumber(start1,start2));
				break;
			case genpafShingle:
				diag = diagNumber (start1, start2);
				diagSE = seq1Len - diag;
				diagNW = seq2Len + diag;
				if (diag < 0)
					{
					if ((diagNW < 0)
					 || ((u32) diagNW < seq1Len)) diag = -diagNW;
					                         else diag = 0;
					}
				else if (diag > 0)
					{
					if ((diagSE < 0)
					 || ((u32) diagSE < seq2Len)) diag = diagSE;
					                         else diag = 0;
					}
				if (diag == 0) fprintf (f, "NA");
				          else fprintf (f, sgnposFmt, diag);
				break;
			case genpafScore:
				fprintf (f, scoreFmt, s);
				break;
			case genpafIdentity:
				segment_identity (seq1, pos1, seq2, pos2, length, &numer, &denom);
				fprintf (f, unsposSlashFmt, numer, denom);
				if (denom != 0) fprintf (f, "\t%.1f%%", (100.0*numer) / denom);
				           else fprintf (f, "\tNA");
				break;
			case genpafCoverage:
				seg.pos1   = pos1;
				seg.pos2   = pos2;
				seg.length = length;
				segment_coverage (seq1, seq2, &seg, &numer, &denom);
				fprintf (f, unsposSlashFmt, numer, denom);
				if (denom != 0) fprintf (f, "\t%.1f%%", (100.0*numer) / denom);
				           else fprintf (f, "\tNA");
				break;
			case genpafContinuity:
				numer = denom = length;
				fprintf (f, unsposSlashFmt, numer, denom);
				if (denom != 0) fprintf (f, "\t%.1f%%", (100.0*numer) / denom);
				           else fprintf (f, "\tNA");
				break;
			case genpafGapRate:
				numer = 0;
				denom = length;
				fprintf (f, unsposSlashFmt, numer, denom);
				if (denom != 0) fprintf (f, "\t%.1f%%", (100.0*numer) / denom);
				           else fprintf (f, "\tNA");
				break;
			default:
				break;
			}
		}
	fprintf (f, "\n");
	}

//----------
//
// parse_genpaf_keys--
//	Convert a list of genpaf names to field keys.
//
//----------
//
// Arguments:
//	char*	s:	A comma-separated list of the NAMES of the fields to print
//	            (genpafName or genpafAliases strings).
//
// Returns:
//	A list of the fields to print (genpafXXX values).  This is dynamically
//	allocated from the heap.
//
//----------

char* parse_genpaf_keys
   (char*	s)
	{
	char*	ss, *keys, *kScan, *field, *diffChars;
	int		numFields, ix, key;
	char	terminator;
	int		haveDiff;

	// figure out how many fields there will be

	numFields = 1;
	for (ss=s ; *ss!=0 ; ss++)
		if (*ss == ',') numFields++;

	// allocate key string (the extra characters are for a potential set of
	// characters for genpafTextDiff)

	keys = malloc_or_die ("parse_genpaf_keys",
	                      numFields+genpafTDInfoSize+1);
	diffChars = NULL;
	haveDiff  = false;

	// parse fields

	kScan = keys;
	field = s;
	for (ss=s ; ; ss++)
		{
		if ((*ss != ',') && (*ss != 0)) continue;

		terminator = *ss;
		*ss = 0;

		if (field[0] == 0) // (empty string)
			key = genpafCR;
		else
			{
			key = -1;
			for (ix=0 ; genpafName[ix].name!=NULL ; ix++)
				{
				if (strcmp (field, genpafName[ix].name) != 0) continue;
				key = genpafName[ix].key;
				break;
				}
			if (key < 0)
				{
				for (ix=0 ; genpafAliases[ix].name!=NULL ; ix++)
					{
					if (strcmp (field, genpafAliases[ix].name) != 0) continue;
					key = genpafAliases[ix].key;
					break;
					}
				}
			if ((key < 0)
			 && (strcmp_prefix (field, genpafTDName) == 0)
			 && (strlen (field) == strlen(genpafTDName) + genpafTDInfoSize))
				{
				key = genpafTextDiff;
				diffChars = field + strlen(genpafTDName);
				if (strchr (diffChars, genpafInfoSeparator) != NULL)
					suicidef ("%s field cannot contain %c",
							  genpafTDName, genpafInfoSeparator);
				}
			if (key < 0)
				suicidef ("unrecognized field name: %s", field);
			if (key == genpafTextDiff)
				{
				if (haveDiff)
					suicidef ("duplicate field name: %s", genpafTDName);
				haveDiff = true;
				}
			}

		*(kScan++) = key;
		field = ss+1;

		if (terminator == 0) break;
		}

	if (diffChars != NULL)
		{
		*(kScan++) = genpafInfoSeparator;
		*(kScan++) = genpafTextDiff;
		for (ix=0 ; ix<genpafTDInfoSize ; ix++)
			*(kScan++) = *(diffChars++);
		}

	*kScan = 0;

	return keys;
	}

//----------
//
// extract_key_info--
//	Find info for a specific field in a list of genpaf field keys.
//
//----------
//
// Arguments:
//	char*	keys:		A list of the fields (genpafXXX values).
//	char	desiredKey:	The key desired (a genpafXXX value).
//
// Returns:
//	A pointer to the key's info, if any;  NULL otherwise.  Note that the pointer
//	is into the keys string, and the info string may be followed by info for
//	other keys.
//
//----------

static char* extract_key_info
   (char*	keys,
	char	desiredKey)
	{
	char*	k;

	for (k=keys ; (*k!=0)&&(*k!=genpafInfoSeparator); k++)
		;

	if (*k == 0) return NULL;

	for ( ; *k!=0; k++)
		{
		k++;								// (step past genpafInfoSeparator) 
		if (*k == desiredKey) return k+1;	// (next character is the key)
		k++;								// (step past the key)
		for ( ; (*k!=0)&&(*k!=genpafInfoSeparator); k++)
			;								// (skip to next genpafInfoSeparator)
		}

	return false;
	}

