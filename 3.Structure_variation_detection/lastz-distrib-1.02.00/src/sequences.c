//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: sequences.c
//
//----------
//
// sequences--
//	Support for DNA sequences.
//
// FASTA format stores DNA sequences as plain text.  A single file can store
// multiple sequences, which we call contigs.  Each contig begins with a header
// line, which begins with a ">".
//
// NIB format stores a single DNA sequence, containing {A,C,G,T,a,c,g,t,N,N} in
// two bases per byte.  As of Jan/2008, a spec for NIB files can be found at
//	http://genome.ucsc.edu/FAQ/FAQformat#format8
//
// 2BIT format stores multiple DNA sequences encoded as four bases per byte
// with some additional information describing runs of masked bases or Ns.  As
// of Jan/2008, a spec for 2BIT files can be found at
//	http://genome.ucsc.edu/FAQ/FAQformat#format7
//
// HSX format is a hashed sequence index that consists of references to
// sequences in other files.  This allows random access.  See the file format
// spec in the lastz readme file for more information.
//
// Quantum-dna format stores each base as a byte, with the meaning of the
// byte values essentially defined by the scoring matrix.  See the file format
// spec in the lastz readme file for more information.
//
//----------
//
// WARNING:	As of this writing (Apr/2008), the code to read sequences is a
//			mess.  The additions of rewindability and contigs-of-interest did
//			not fit well with the original design, and have been patched in
//			with bandaids.  The author hopes to rectify this in the future.
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
#include <ctype.h>				// standard C upper/lower stuff
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff

#define  sequences_owner		// (make this the owner of its globals)
#include "sequences.h"			// interface to this module


#define pathSlash '/'
#ifdef compileForWindows
#undef pathSlash
#define pathSlash '\\'
#endif

//----------
//
// debug set ups
//
//----------

#define sequence_filename(seq)  (((seq)->fileName != NULL)? (seq)->fileName : "(unnamed sequence file)")

//--- debug set up for debugging the contigs-of-interest names file ---

//#define debugNamesFile

#ifndef debugNamesFile
#define debugNamesFile_1 ;
#define debugNamesFile_2 ;
#define debugNamesFile_3 ;
#define debugNamesFile_4 ;
#define debugNamesFile_5 ;
#define debugNamesFile_6 ;
#define debugNamesFile_7 ;
#define debugNamesFile_8 ;
#endif // not debugNamesFile

#ifdef debugNamesFile

#define debugNamesFile_1                                                     \
	fprintf (stderr, "load_sequence (%s):\n", sequence_filename(_seq));      \
	if ((_seq != NULL) && (_seq->preLoaded))                                 \
		fprintf (stderr, "  preloaded, header: %s\n", _seq->header);

#define debugNamesFile_2                                                     \
	fprintf (stderr, "  header: %s\n", _seq->header);

#define debugNamesFile_3                                                     \
	fprintf (stderr, "another_sequence (%s):\n", sequence_filename(_seq));

#define debugNamesFile_4                                                     \
	long int	fpos;                                                        \
	fprintf (stderr, "find_next_general_fasta_coi (%s):\n", sequence_filename(_seq));

#define debugNamesFile_5                                                     \
		fpos = ftell (_seq->f) - _seq->pendingLen - 1;

#define debugNamesFile_6                                                     \
	fprintf (stderr, "  found: [%08lX] %s\n", fpos, _seq->nextContigName);

#define debugNamesFile_7                                                     \
	fprintf (stderr, "find_next_2bit_coi (%s): [%s]\n",                      \
	                 sequence_filename(_seq), _seq->nextContigName);

#define debugNamesFile_8                                                     \
	fprintf (stderr, "find_next_hsx_coi (%s): [%s]\n",                       \
	                 sequence_filename(_seq), _seq->nextContigName);

#endif // debugNamesFile

//--- debug set up for debugging partitioned sequences ---

//#define debugPartitions

#ifndef debugPartitions
#define debugPartitions_1 ;
#define debugPartitions_2 ;
#endif // not debugPartition

#ifdef debugPartitions

static void print_partition_table (FILE* f, seq* _seq);

#define debugPartitions_1                                                   \
	_seq->v[_seq->len] = '*';

#define debugPartitions_2                                                   \
	if (sp->p != NULL)                                                      \
		{                                                                   \
		print_partition_table (stderr, _seq);                               \
		print_sequence (stderr, _seq, "", 100);                             \
		}

#endif // debugPartitions

//----------
//
// private global data relating to fasta and csfasta format
//
//----------

// tables to map 8-bit ascii character to fasta or csfasta character type
//   "nucleotide" characters are the A, C, G, T and N
//   "ambiguous"  characters are the remaining IUPAC 15-letter alphabet

enum
	{ _bad = 0, _whitespace, _newline, _nucleotide, _ambiguous, _color };

#define __ _bad
#define _w _whitespace
#define _l _newline
#define _n _nucleotide
#define _a _ambiguous
#define _c _color

static const u8 char_to_fasta_type[256] =
	{
	__,__,__,__,__,__,__,__,__,_w,_l,__,_w,_l,__,__, // 0x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 1x
	_w,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 2x
	_w,_w,_w,_w,_w,_w,_w,_w,_w,_w,__,__,__,__,__,__, // 3x (numbers)
	__,_n,_a,_n,_a,__,__,_n,_a,__,__,_a,__,_a,_n,__, // 4x (upper case)
	__,__,_a,_a,_n,__,_a,_a,_n,_a,__,__,__,__,__,__, // 5x (upper case)
	__,_n,_a,_n,_a,__,__,_n,_a,__,__,_a,__,_a,_n,__, // 6x (lower case)
	__,__,_a,_a,_n,__,_a,_a,_n,_a,__,__,__,__,__,__, // 7x (lower case)
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 8x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 9x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ax
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Bx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Cx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Dx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ex
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__  // Fx
	};

static const u8 char_to_csfasta_type[256] =
	{
	__,__,__,__,__,__,__,__,__,_w,_l,__,_w,_l,__,__, // 0x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 1x
	_w,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 2x
	_c,_c,_c,_c,__,__,__,__,__,__,__,__,__,__,__,__, // 3x (numbers)
	__,_n,__,_n,__,__,__,_n,__,__,__,__,__,__,__,__, // 4x (upper case)
	__,__,__,__,_n,__,__,__,__,__,__,__,__,__,__,__, // 5x (upper case)
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 6x (lower case)
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 7x (lower case)
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 8x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 9x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ax
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Bx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Cx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Dx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ex
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__  // Fx
	};

#undef __
#undef _w
#undef _l
#undef _n
#undef _a
#undef _c

//----------
//
// private global data relating to nib format
//
//----------

// nib file magic number(s)

static const u32 nibMagicBig     = 0x6BE93D3A;	// in big endian format
static const u32 nibMagicLittle  = 0x3A3DE96B;	// in little endian format

// tables to map 4-bit nybbles from nib file to a character
//
// nybbles map as follows:
//		nybble:    0 1 2 3 4 5 6 7 8 9 A B C D E F
//		character: T C A G N X X X t c a g n x x x
// For (alleged) efficiency's sake, we use separate lookup tables for the left
// and right nybble mapping.

static const unsigned char nibTo1stChar[256] = 
	"TTTTTTTTTTTTTTTT"
	"CCCCCCCCCCCCCCCC"
	"AAAAAAAAAAAAAAAA"
	"GGGGGGGGGGGGGGGG"
	"NNNNNNNNNNNNNNNN"
	"XXXXXXXXXXXXXXXX"
	"XXXXXXXXXXXXXXXX"
	"XXXXXXXXXXXXXXXX"
	"tttttttttttttttt"
	"cccccccccccccccc"
	"aaaaaaaaaaaaaaaa"
	"gggggggggggggggg"
	"nnnnnnnnnnnnnnnn"
	"xxxxxxxxxxxxxxxx"
	"xxxxxxxxxxxxxxxx"
	"xxxxxxxxxxxxxxxx";

static const unsigned char nibTo2ndChar[256] = 
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx";

static const unsigned char nibTo1stCharUnmasked[256] = 
	"TTTTTTTTTTTTTTTT"
	"CCCCCCCCCCCCCCCC"
	"AAAAAAAAAAAAAAAA"
	"GGGGGGGGGGGGGGGG"
	"NNNNNNNNNNNNNNNN"
	"XXXXXXXXXXXXXXXX"
	"XXXXXXXXXXXXXXXX"
	"XXXXXXXXXXXXXXXX"
	"TTTTTTTTTTTTTTTT"
	"CCCCCCCCCCCCCCCC"
	"AAAAAAAAAAAAAAAA"
	"GGGGGGGGGGGGGGGG"
	"NNNNNNNNNNNNNNNN"
	"XXXXXXXXXXXXXXXX"
	"XXXXXXXXXXXXXXXX"
	"XXXXXXXXXXXXXXXX";

static const unsigned char nibTo2ndCharUnmasked[256] = 
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX";

//----------
//
// private global data relating to 2bit format
//
//----------

// 2bit file magic number(s)

static const u32 twobitMagicBig    = 0x1A412743;	// in big endian format
static const u32 twobitMagicLittle = 0x4327411A;	// in little endian format

static const char* twobitToChars[256] =
	{
	"TTTT","TTTC","TTTA","TTTG","TTCT","TTCC","TTCA","TTCG",
	"TTAT","TTAC","TTAA","TTAG","TTGT","TTGC","TTGA","TTGG",
	"TCTT","TCTC","TCTA","TCTG","TCCT","TCCC","TCCA","TCCG",
	"TCAT","TCAC","TCAA","TCAG","TCGT","TCGC","TCGA","TCGG",
	"TATT","TATC","TATA","TATG","TACT","TACC","TACA","TACG",
	"TAAT","TAAC","TAAA","TAAG","TAGT","TAGC","TAGA","TAGG",
	"TGTT","TGTC","TGTA","TGTG","TGCT","TGCC","TGCA","TGCG",
	"TGAT","TGAC","TGAA","TGAG","TGGT","TGGC","TGGA","TGGG",

	"CTTT","CTTC","CTTA","CTTG","CTCT","CTCC","CTCA","CTCG",
	"CTAT","CTAC","CTAA","CTAG","CTGT","CTGC","CTGA","CTGG",
	"CCTT","CCTC","CCTA","CCTG","CCCT","CCCC","CCCA","CCCG",
	"CCAT","CCAC","CCAA","CCAG","CCGT","CCGC","CCGA","CCGG",
	"CATT","CATC","CATA","CATG","CACT","CACC","CACA","CACG",
	"CAAT","CAAC","CAAA","CAAG","CAGT","CAGC","CAGA","CAGG",
	"CGTT","CGTC","CGTA","CGTG","CGCT","CGCC","CGCA","CGCG",
	"CGAT","CGAC","CGAA","CGAG","CGGT","CGGC","CGGA","CGGG",

	"ATTT","ATTC","ATTA","ATTG","ATCT","ATCC","ATCA","ATCG",
	"ATAT","ATAC","ATAA","ATAG","ATGT","ATGC","ATGA","ATGG",
	"ACTT","ACTC","ACTA","ACTG","ACCT","ACCC","ACCA","ACCG",
	"ACAT","ACAC","ACAA","ACAG","ACGT","ACGC","ACGA","ACGG",
	"AATT","AATC","AATA","AATG","AACT","AACC","AACA","AACG",
	"AAAT","AAAC","AAAA","AAAG","AAGT","AAGC","AAGA","AAGG",
	"AGTT","AGTC","AGTA","AGTG","AGCT","AGCC","AGCA","AGCG",
	"AGAT","AGAC","AGAA","AGAG","AGGT","AGGC","AGGA","AGGG",

	"GTTT","GTTC","GTTA","GTTG","GTCT","GTCC","GTCA","GTCG",
	"GTAT","GTAC","GTAA","GTAG","GTGT","GTGC","GTGA","GTGG",
	"GCTT","GCTC","GCTA","GCTG","GCCT","GCCC","GCCA","GCCG",
	"GCAT","GCAC","GCAA","GCAG","GCGT","GCGC","GCGA","GCGG",
	"GATT","GATC","GATA","GATG","GACT","GACC","GACA","GACG",
	"GAAT","GAAC","GAAA","GAAG","GAGT","GAGC","GAGA","GAGG",
	"GGTT","GGTC","GGTA","GGTG","GGCT","GGCC","GGCA","GGCG",
	"GGAT","GGAC","GGAA","GGAG","GGGT","GGGC","GGGA","GGGG"
	};

//----------
//
// private global data relating to hsx format
//
//----------

// hsx file magic number(s)

static const u32 hsxMagicBig    = 0xD2527095;	// in big endian format
static const u32 hsxMagicLittle = 0x957052D2L;	// in little endian format

static const u64 hsxMsBit5      = ((u64) 0x80) << (4*8);
static const u64 hsxMaxFilePos  = (u64) ((long int) -1);

//----------
//
// private global data relating to quantum-dna format
//
//----------

// quantum-dna file magic number(s)

static const u32 qdnaMagicBig       = 0xC4B47197;	// in big endian format
static const u32 qdnaMagicLittle    = 0x9771B4C4;	// in little endian format
static const u32 oldQdnaMagicBig    = 0xF656659E;	// in big endian format
static const u32 oldQdnaMagicLittle = 0x9E6556F6;	// in little endian format

//----------
//
// prototypes for private functions
//
//----------

static seq*   alloc_sequence_record (char* id);
static void   skip_sequences        (seq* _seq, int skipCount);
static void   load_sequence_core    (seq* _seq, int keeper);
static void   load_fasta_sequence   (seq* _seq, int keeper);
static unspos parse_fasta           (seq* _seq, int storeEm);
static void   load_csfasta_sequence (seq* _seq, int keeper);
static unspos parse_csfasta         (seq* _seq, int storeEm);
static void   load_nib_sequence     (seq* _seq, int keeper);
static void   read_2bit_header      (seq* _seq);
static void   load_2bit_sequence    (seq* _seq, int keeper);
static void   read_hsx_header       (seq* _seq);
static void   load_hsx_sequence     (seq* _seq, int keeper);
static void   load_qdna_sequence    (seq* _seq, int keeper);
static int    another_sequence_core (seq* _seq);
static void   create_short_header   (seq* _seq);
static int    find_next_fasta_coi   (seq* _seq);
static int    find_next_csfasta_coi (seq* _seq);
static int    find_next_general_fasta_coi (seq* _seq, int allowComments);
static int    find_next_2bit_coi    (seq* _seq);
static int    find_next_hsx_coi     (seq* _seq);
static int    read_contig_name      (seq* _seq);
static void   shorten_header        (char* src, int nameParseType, int skipPath,
                                     char** dst, u32* dstSize);
static void   expand_nickname       (char* src, u32 contigNumber,
                                     char** dst, u32* dstSize);
static void   enough_partitions     (seq* _seq, u32 numPartitions, u32 poolSize,
                                     int anticipate);
static void   parse_sequence_name   (const char* name,
                                     char** fileName, char** nickname,
                                     char** contigOfInterest,
                                     char** namesFileName,
                                     int* subsampleK, int* subsampleN,
                                     char** softMaskFileName, int* softMaskComplement,
                                     char** xMaskFileName, int* xMaskComplement,
                                     char** nMaskFileName, int* nMaskComplement,
                                     int* nameParseType,
                                     char** nameTrigger,
                                     int* doRevCompFlags,
                                     int* doUnmask, int* doJoin,
                                     int* useFullNames,
                                     int* fileType,
                                     int* isQuantum, char** qCodingFileName,
                                     unspos* _start, unspos* _end,
                                     int* endIsSoft);
static int    detect_file_type      (seq* _seq);
static u32    read_4                (seq* _seq, int asBigEndian);
static u32    read_4_big            (seq* _seq);
static u32    read_4_little         (seq* _seq);
static u64    read_5                (seq* _seq, int asBigEndian);
static u64    read_5_big            (seq* _seq);
static u64    read_5_little         (seq* _seq);
static u64    read_6                (seq* _seq, int asBigEndian);
static u64    read_6_big            (seq* _seq);
static u64    read_6_little         (seq* _seq);
static int    skip_seq_whitespace   (seq* _seq);
static int    skip_fasta_whitespace (seq* _seq);
static int    seq_getc              (seq* _seq);
static void   seq_ungetc            (char ch, seq* _seq);
static int    skip_chars            (seq* _seq, u32 toSkip);
static int    test_rewindability    (seq* _seq);
static void   save_fstate           (seq* _seq);
static void   restore_fstate        (seq* _seq);
static void   dump_sequence         (FILE* f, seq* _seq);

//----------
//
// open_sequence_file--
//	Open a sequence file for read operations.
// open_rewindable_sequence_file--
//	Open a sequence file for read operations, which may be rewound later.  Note
//	that if the actual file is not rewindable, but only contains one sequence,
//	we can still satisfy the caller's desire for rewindability.
//
//----------
//
// Arguments:
//	char*	name:			The name of the file that sequence data is to be
//							.. read from.  This can be NULL, which indicates
//							.. that data is to be read from stdin.
// 	int		fileType:		The type of file, e.g. seq_type_nib.  This can be
//							.. seq_type_unknown if the caller would like the
//							.. type to be determined from the file.
//	int		needTrueLen:	true  => set seq->trueLen correctly, even if this
//							         .. means reading additional characters
//							         .. outside the desired (sub)interval
//							false => the value written to trueLen is unimportant
//	int		allowAmbiDNA:	true  => permit ambiguous DNA characters
//							          .. B,D,H,K,M,R,S,V,W,Y
//							false => only A,C,G,T,N,X permitted
//	u8*		qToComplement:	(similar to nuc_to_complement) array to map a
//							.. quantum base to its complement.  This is only
//							.. used if the sequence is quantum DNA.  This may
//							.. be NULL (in which case the sequence can not be
//							.. reverse-complemented).
//
// Returns:
//	A pointer to the sequence;  failures result in fatality.  The caller must
//	eventually de-allocate this by calling free_sequence().
//
//----------
//
// Implementation notes:
//
// To satisfy the caller's request that the file be rewindable, we perform the
// initial sequence load here (and set a flag to let load_sequence know this
// has happened).  Then we check whether the file contains any additional data.
// If it doesn't, then we treat the file as rewindable regardless of whether
// the underlying file actualy is.  Only if the file contains additional data
// do we then perform a test for rewindability, by attempting to set the file's
// position back to the end of the first sequence.
//
//----------

static seq* private_open_sequence_file (char* name, int fileType,
               int beRewindable, int needTrueLen, int allowAmbiDNA,
               u8* qToComplement);

seq* open_sequence_file (char* name, int fileType, int needTrueLen, int allowAmbiDNA, u8* qToComplement)
	{ return private_open_sequence_file (name, fileType, false, needTrueLen, allowAmbiDNA, qToComplement); }

seq* open_rewindable_sequence_file (char* name, int fileType, int needTrueLen, int allowAmbiDNA, u8* qToComplement)
	{ return private_open_sequence_file (name, fileType, true, needTrueLen, allowAmbiDNA, qToComplement); }

static seq* private_open_sequence_file
   (char*	name,
	int		fileType,
	int		beRewindable,
	int		needTrueLen,
	int		allowAmbiDNA,
	u8*		qToComplement)
	{
	seq*	_seq;
	int		isQuantum = false;
	char*	qCodingFileName = NULL;
	int		forcedfileType = seq_type_unknown;
	int		err;

	// allocate the sequence tracking structure

	_seq = alloc_sequence_record ("open_sequence");
	_seq->vOwner  = true;  				// (even though _seq->v is NULL, we
	_seq->vcOwner = true;				//  .. still will be considered as the
										//  .. 'owner' so we can resize it)
	_seq->partition.poolOwner = true;	// (similarly, we are owner of partition
										//  .. names pool so we can resize it)
	_seq->headerOwner         = true;  	// (and we're owner of header[] and
	_seq->shortHeaderOwner    = true;	//  .. shortHeader[] so we can resize
										//  .. them)

	_seq->pendingChars = zalloc_or_die ("open_sequence (pendingChars)",
                                        seqBufferSize);
	_seq->pendingStack = _seq->pendingChars + seqBufferSize;
	_seq->pendingLen   = 0;

	_seq->lockedHeader = false;
	_seq->needTrueLen  = needTrueLen;
	_seq->allowAmbiDNA = allowAmbiDNA;

	// if we have no name, use stdin

	if (name == NULL)
		{
		_seq->fileName = copy_string ("(stdin)");
		_seq->f        = stdin;
		}

	// otherwise, open the file
	// nota bene:  we'd like to open the file as "rt" instead of "rb" if it is
	//             a fasta file;  unfortunately we don't know what the file
	//             type is until we open it

	else
		{
		parse_sequence_name (name,
		                     &_seq->fileName, &_seq->header,
		                     &_seq->contigOfInterest,
		                     &_seq->namesFileName,
                             &_seq->subsampleK, &_seq->subsampleN,
		                     &_seq->softMaskFileName, &_seq->softMaskComplement,
		                     &_seq->xMaskFileName, &_seq->xMaskComplement,
		                     &_seq->nMaskFileName, &_seq->nMaskComplement,
		                     &_seq->nameParseType,
		                     &_seq->nameTrigger,
		                     &_seq->doRevCompFlags,
		                     &_seq->doUnmask, &_seq->doJoin,
		                     &_seq->useFullNames,
		                     &forcedfileType,
		                     &isQuantum, &qCodingFileName,
		                     &_seq->startLimit, &_seq->endLimit,
		                     &_seq->endIsSoft);
		_seq->f = fopen_or_die (_seq->fileName, "rb");

		if (_seq->header != NULL)
			{
			_seq->headerSize   = strlen (_seq->header) + 1;
			_seq->lockedHeader = true;
			_seq->hasNickname  = true;
			}
		}

	// init any non-zero fields

	_seq->hasSavedState = false;
	_seq->rewindable    = -1;	// (rewindability unknown at this point)
	_seq->contig        = 0;

	if (isQuantum)
		{
		if ((fileType != seq_type_qdna) && (fileType != seq_type_unknown))
			suicidef ("clashing file type for %s ([quantum] used for %s file)",
			          sequence_filename(_seq), seqTypeNames[fileType]);
		_seq->fileType = seq_type_qdna;
		}
	else if (forcedfileType != seq_type_unknown)
		{
		if ((fileType != forcedfileType) && (fileType != seq_type_unknown))
			suicidef ("clashing file type for %s (%s used for %s file)",
			          sequence_filename(_seq), seqTypeNames[forcedfileType], seqTypeNames[fileType]);
		_seq->fileType = forcedfileType;
		}
	else
		_seq->fileType = fileType;

	// initialize subsampling

	if (_seq->subsampleN == 1)
		_seq->subsampleN = 0; // (subsampling 1 of 1 is meaningless)

	if (_seq->subsampleN == 0)
		_seq->subsampleSkip = 0;
	else
		_seq->subsampleSkip = _seq->subsampleK - 1;

	// if the sequences in this file are to be joined into a parititioned
	// sequence, initialize that

	if (_seq->doJoin)
		{
		enough_partitions (_seq, /*numPartitions*/ 100, /*poolSize*/ 0,
		                   /*anticipate*/ false);
		_seq->partition.state = seqpart_empty;
		}

	// if the file type is not yet known, figure out what it is

	if (_seq->fileType == seq_type_unknown)
		_seq->fileType = detect_file_type (_seq);

	if ((!sequences_dbgAllowColors)
	 && (_seq->fileType == seq_type_csfasta))
		suicidef ("sorry, color space is not fully implemented yet");

	if ((_seq->fileType != seq_type_2bit)
	 && (_seq->fileType != seq_type_hsx)
	 && (_seq->contigOfInterest != NULL))
		suicidef ("specific contig-of-interest only valid for 2bit or hsx files (%s)",
		          _seq->contigOfInterest);

	if ((_seq->fileType != seq_type_fasta)
	 && (_seq->fileType != seq_type_csfasta)
	 && (_seq->fileType != seq_type_2bit)
	 && (_seq->fileType != seq_type_hsx)
	 && (_seq->namesFileName != NULL))
		suicidef ("sequence-subset file only valid for fasta, csfasta, 2bit or hsx files\n(%s)",
		          _seq->namesFileName);

	if ((_seq->fileType == seq_type_hsx)
	 && (_seq->nameParseType != name_parse_type_core))
		suicidef ("\"nameparse=\" is not valid for hsx files");

	if ((_seq->fileType != seq_type_fasta)
	 && (_seq->fileType != seq_type_csfasta)
	 && (_seq->fileType != seq_type_2bit)
	 && (_seq->nameParseType == name_parse_type_alnum))
		suicidef ("\"nameparse=alphanum\" only valid for fasta, csfasta or 2bit files");

	if ((_seq->fileType != seq_type_fasta)
	 && (_seq->fileType != seq_type_csfasta)
	 && (_seq->fileType != seq_type_2bit)
	 && (_seq->nameParseType == name_parse_type_darkspace))
		suicidef ("\"nameparse=darkspace\" only valid for fasta, csfasta or 2bit files");

	if ((_seq->fileType != seq_type_fasta)
	 && (_seq->fileType != seq_type_csfasta)
	 && (_seq->nameTrigger != NULL))
		suicidef ("\"nameparse=tag:%s\" only valid for fasta or csfasta files", _seq->nameTrigger);

	// for quantum DNA files, attach the complement mapping

	if (_seq->fileType != seq_type_qdna)
		_seq->qToComplement = qToComplement;

	// make sure the file is rewindable if the caller requires it to be

	if (beRewindable)
		{
		if (load_sequence (_seq))
			{
			_seq->preLoaded = true;
			if (another_sequence (_seq))
				{
				err = test_rewindability (_seq);
				if (err != 0) goto not_rewindable;
				_seq->rewindable = true;
				}
			}
		}

	// for 2bit or hsx files, we need to read the header
	// $$$ should this be moved to before the pre-load above?

	if (_seq->fileType == seq_type_2bit)
		read_2bit_header (_seq);
	else if (_seq->fileType == seq_type_hsx)
		read_hsx_header (_seq);

	// if we have a contigs-of-interest file, open it, read the first line, and
	// advance to that contig

	_seq->contigPending = false;

	if (_seq->namesFileName != NULL)
		{
		_seq->namesFile = fopen_or_die (_seq->namesFileName, "rt");
		if (!read_contig_name (_seq))
			suicidef ("contigs-of-interest file is empty: %s", _seq->namesFile);

		if (_seq->fileType == seq_type_fasta)
			find_next_fasta_coi (_seq);
		else if (_seq->fileType == seq_type_csfasta)
			find_next_csfasta_coi (_seq);
		else if (_seq->fileType == seq_type_2bit)
			find_next_2bit_coi (_seq);
		else // if (_seq->fileType == seq_type_hsx)
			find_next_hsx_coi (_seq);
		}

	// if we have a quantum coding file, read it

	if (qCodingFileName != NULL)
		{
		_seq->qCoding = read_quantum_code_by_name (qCodingFileName);
		free_if_valid   ("open_sequence_file (qCodingFileName)", qCodingFileName);
		}

	return _seq;

// failure exits

not_rewindable:
	if (name == NULL) name = "(stdin)";
	suicidef_with_perror ("sequence file \"%s\" is not rewindable"
						  " (fseek returned %d): %s",
						  name, err, sequence_filename(_seq));
	return NULL; // (never gets here)
	}

//----------
//
// rewind_sequence_file--
//	Rewind a sequence file, so that the sequence(s) within it can be read again.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to rewind.
//
// Returns:
//	(nothing;  failures cause program termination)
//
//----------

void rewind_sequence_file
   (seq*	_seq)
	{
	if (_seq->f == NULL)  return;	// (no file to rewind)

	// decide whether we can take a short cut and ignore the rewind request

	if ((_seq->fileType == seq_type_2bit)		// (we have only one
	 && (_seq->contigOfInterest != NULL))		//  .. sequence in the file)
		{
		if (_seq->contig != 0) _seq->preLoaded = true;
		return;
		}

	if ((_seq->fileType == seq_type_hsx)		// (we have only one
	 && (_seq->contigOfInterest != NULL))		//  .. sequence in the file)
		{
		if (_seq->contig != 0) _seq->preLoaded = true;
		return;
		}

	if (_seq->contig < 2)						// (we haven't read more than
		{										//  .. one sequence from file)
		if (_seq->contig != 0) _seq->preLoaded = true;
		return;
		}

	if (_seq->fileType == seq_type_2bit)
		{
		_seq->twoBit.contigFilePos = _seq->twoBit.indexFilePos;
		_seq->twoBit.contigLoaded  = false;
		goto reset_file_data;
		}

	if (_seq->fileType == seq_type_hsx)
		{
		_seq->hsx.contigLoaded = false;
		goto reset_file_data;
		}

	if (_seq->rewindable == -1)
		_seq->rewindable = (test_rewindability (_seq) == 0);

	if (_seq->rewindable == false)
		suicidef ("sequence file is not rewindable: %s", sequence_filename(_seq));

	// rewind the file and reset the data

	rewind (_seq->f);

reset_file_data:
	_seq->len           = 0;
	_seq->contig        = 0;
	_seq->preLoaded     = false;
	_seq->pendingStack  = _seq->pendingChars + seqBufferSize;
	_seq->pendingLen    = 0;
	_seq->hasSavedState = false;

	// rewind the contigs-of-interest file

	if (_seq->namesFile != NULL)
		{
		int err = fseek (_seq->namesFile, 0, SEEK_SET);
		if (err != 0)
			suicidef ("failed to seek to position in file\n"
			          "in rewind_sequence_file for %s, index fseek(%08lX) returned %d",
			          _seq->namesFileName, 0, err);

		read_contig_name (_seq);
		}

	}

//----------
//
// copy_sequence--
//	Make a copy of a sequence.  Note that the copy is not as functional as the
//	original.  For example, file operations cannot be performed on the copy.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to make a copy of.
//
// Returns:
//	A pointer to the sequence;  failures result in fatality.  The caller must
//	eventually de-allocate this by calling free_sequence().
//
//----------

seq* copy_sequence
   (seq*	_seq)
	{
	seq*	newSeq;
	u32		ix;

	// allocate the sequence tracking structure

	newSeq = alloc_sequence_record ("copy_sequence");

	newSeq->pendingChars = zalloc_or_die ("copy_sequence (pendingChars)",
	                                      seqBufferSize);
	newSeq->pendingStack = newSeq->pendingChars + seqBufferSize;
	newSeq->pendingLen   = 0;

	// allocate the sequence content

	if ((_seq->v == NULL) || (_seq->len < 1))
		{
		newSeq->v       = NULL;
		newSeq->vOwner  = true;
		newSeq->size    = 0;
		newSeq->len     = 0;
		newSeq->trueLen = 0;
		}
	else
		{
		if (_seq->len > maxSequenceLen)
			suicidef ("in copy_sequence, "
			          "sequence length " unsposFmt " exceeds maximum (" unsposFmt ")",
			          _seq->len, maxSequenceLen);

		newSeq->v       = malloc_or_die ("copy_sequence (v)", _seq->len+1);
		newSeq->vOwner  = true;
		newSeq->size    = _seq->len+1;
		newSeq->len     = _seq->len;
		newSeq->trueLen = _seq->trueLen;
		}

	if (_seq->vc == NULL)
		{
		newSeq->vc      = NULL;
		newSeq->vcOwner = true;
		}
	else
		{
		if (_seq->len > maxSequenceLen)
			suicidef ("in copy_sequence, "
			          "sequence length " unsposFmt " exceeds maximum (" unsposFmt ")",
			          _seq->len, maxSequenceLen);

		newSeq->vc      = malloc_or_die ("copy_sequence (v)", _seq->len+1);
		newSeq->vcOwner = true;
		newSeq->len     = _seq->len;
		}

	// set up file info

	newSeq->fileType = seq_type_nofile;

	if (_seq->fileName != NULL)
		newSeq->fileName = copy_string (_seq->fileName);

	// copy the sequence content (including a terminating zero)

	if (newSeq->v != NULL)
		{
		for (ix=0 ; ix<=newSeq->len ; ix++)
			newSeq->v[ix] = _seq->v[ix];
		}

	newSeq->preLoaded = _seq->preLoaded;
	newSeq->contig    = _seq->contig;

	// copy other fields

	newSeq->start        = _seq->start;
	newSeq->revCompFlags = _seq->revCompFlags;
	newSeq->contig       = _seq->contig;
	newSeq->lockedHeader = _seq->lockedHeader;
	newSeq->allowAmbiDNA = _seq->allowAmbiDNA;

	if (_seq->header == NULL)
		{
		newSeq->header     = NULL;
		newSeq->headerSize = 0;
		}
	else
		{
		newSeq->header     = copy_string (_seq->header);
		newSeq->headerSize = strlen (newSeq->header) + 1;
		}
	newSeq->headerOwner = true;

	if (_seq->shortHeader == NULL)
		{
		newSeq->shortHeader     = NULL;
		newSeq->shortHeaderSize = 0;
		}
	else
		{
		newSeq->shortHeader      = copy_string (_seq->shortHeader);
		newSeq->shortHeaderSize  = strlen (newSeq->shortHeader) + 1;
		}
	newSeq->shortHeaderOwner = true;

	newSeq->startLimit   = _seq->startLimit;
	newSeq->endLimit     = _seq->endLimit;
	newSeq->endIsSoft    = _seq->endIsSoft;
	newSeq->useFullNames = _seq->useFullNames;

	if (_seq->contigOfInterest == NULL)
		newSeq->contigOfInterest = NULL;
	else
		newSeq->contigOfInterest = copy_string (_seq->contigOfInterest);

	return newSeq;
	}

//----------
//
// new_sequence--
//	Create a new, empty, sequence.
//
//----------
//
// Arguments:
//	unspos	allocLen:	The sequence length to allocate for.  The special
//						.. value seqposInfinity indicates that no memory should
//						.. be allocated for the sequence vector.
//
// Returns:
//	A pointer to the sequence;  failures result in fatality.  The caller must
//	eventually de-allocate this by calling free_sequence().
//
//----------

seq* new_sequence
   (unspos	allocLen)
	{
	seq*	_seq;

	// allocate the sequence tracking structure

	_seq = alloc_sequence_record ("new_sequence");

	_seq->pendingChars = zalloc_or_die ("new_sequence (pendingChars)",
	                                   seqBufferSize);
	_seq->pendingStack = _seq->pendingChars + seqBufferSize;
	_seq->pendingLen   = 0;

	// allocate space for sequence data (but leave it empty)

	if (allocLen == seqposInfinity)
		{
		_seq->v       = NULL;
		_seq->vc      = NULL;
		_seq->vOwner  = false;
		_seq->vcOwner = false;
		_seq->size    = 0;
		_seq->len     = 0;
		}
	else
		{
		_seq->v       = malloc_or_die ("new_sequence (v)", allocLen+1);
		_seq->vc      = NULL;
		_seq->vOwner  = true;
		_seq->vcOwner = true;
		_seq->size    = allocLen+1;
		_seq->len     = 0;
		_seq->v[0]    = 0;
		}

	// initialize the other fields

	_seq->fileType = seq_type_nofile;
	_seq->contig   = 1;
	_seq->start    = 1;
	_seq->trueLen  = 0;

	return _seq;
	}

//----------
//
// alloc_sequence_record--
//	Allocate a new sequence tracking structure, and make sure all pointer
//	fields are NULL.
//
//----------
//
// Arguments:
//	char*	id:	an identifying string to be used when trackMemoryUsage is
//				.. turned on;  this can be NULL.
//
// Returns:
//	A pointer to the sequence;  failures result in fatality.  The caller must
//	eventually de-allocate this by calling free_sequence().
//
//----------

static seq* alloc_sequence_record
   (arg_dont_complain(char* id))
	{
	seq*	_seq;

	_seq = zalloc_or_die (id, sizeof(*_seq));

	_seq->v                   = NULL;
	_seq->pendingChars        = NULL;
	_seq->fileName            = NULL;
	_seq->header              = NULL;
	_seq->shortHeader         = NULL;
	_seq->f                   = NULL;
	_seq->namesFile           = NULL;
	_seq->namesFileName       = NULL;
	_seq->subsampleK          = 0;
	_seq->subsampleN          = 0;
	_seq->softMaskFileName    = NULL;
	_seq->xMaskFileName       = NULL;
	_seq->nMaskFileName       = NULL;
	_seq->nameTrigger         = NULL;
	_seq->contigOfInterest    = NULL;
	_seq->twoBit.nBlockStarts = NULL;
	_seq->twoBit.nBlockSizes  = NULL;
	_seq->twoBit.mBlockstarts = NULL;
	_seq->twoBit.mBlocksizes  = NULL;
	_seq->partition.p         = NULL;
	_seq->partition.pool      = NULL;
	_seq->allowAmbiDNA        = false;
	_seq->qToComplement       = NULL;

	return _seq;
	}

//----------
//
// sequence_long_enough--
//	Make sure a sequence has enough room (including an extra byte for a
//	terminating zero).
//
//----------
//
// Arguments:
//	seq*	_seq:		The sequence to check.
//	unspos	allocLen:	The sequence length to allocate for (not including the
//						.. terminator).
//	int		anticipate:	true  => allocate extra, anticipating the need for more
//						false => don't
//
// Returns:
//	nothing;  the sequence's v[] may be modified;  failures result in fatality.
//
//----------

void sequence_long_enough
   (seq*	_seq,
	unspos	allocLen,
	int		anticipate)
	{
	if (_seq->size >= allocLen+1)
		return;

	allocLen += 2;						// (add space for a terminating zero,
	if (anticipate)						//  .. etc.)
		allocLen += 30 + allocLen / 8;	// anticipatory, grow by about 13%
	allocLen = round_up_16K (allocLen);	// we expect that allocation in
										// .. multiples of 16K is better for
										// .. the heap manager

	if (!_seq->vOwner)
		{
		char* name = (_seq->fileName != NULL)? _seq->fileName
		                                     : _seq->header;
		suicidef ("internal error, attempt to resize external sequence (%s)",
		          name);
		}

	_seq->v    = realloc_or_die ("sequence_long_enough", _seq->v, allocLen);
	_seq->size = allocLen;
	}

//----------
//
// free_sequence--
//	Deallocate a sequence, along with any associated memory or files.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to dispose of.
//
// Returns:
//	(nothing)
//
//----------

void free_sequence
   (seq*	_seq)
	{
	if (_seq == NULL) return;

	if (_seq->vOwner)           free_if_valid ("free_sequence (v)",                _seq->v);
	if (_seq->vcOwner)          free_if_valid ("free_sequence (vc)",               _seq->vc);
	                            free_if_valid ("free_sequence (pendingChars)",     _seq->pendingChars);
	                            free_if_valid ("free_sequence (filename)",         _seq->fileName);
	if (_seq->headerOwner)      free_if_valid ("free_sequence (header)",           _seq->header);
	if (_seq->shortHeaderOwner) free_if_valid ("free_sequence (shortheader)",      _seq->shortHeader);
	                            free_if_valid ("free_sequence (qCoding)",          _seq->qCoding);

	                            free_if_valid ("free_sequence (contigOfInterest)", _seq->contigOfInterest);

	if (_seq->fileType != seq_type_nofile)
		{
		fclose_if_valid (_seq->f);
		fclose_if_valid (_seq->namesFile);
		free_if_valid   ("free_sequence (namesFileName)",     _seq->namesFileName);
		free_if_valid   ("free_sequence (softMaskFileName)",  _seq->softMaskFileName);
		free_if_valid   ("free_sequence (xMaskFileName)",     _seq->xMaskFileName);
		free_if_valid   ("free_sequence (nMaskFileName)",     _seq->nMaskFileName);
		free_if_valid   ("free_sequence (nameTrigger)",       _seq->nameTrigger);
		}

	free_if_valid ("free_sequence (twoBit.nBlockStarts)",     _seq->twoBit.nBlockStarts);
	free_if_valid ("free_sequence (twoBit.nBlockSizes)",      _seq->twoBit.nBlockSizes);
	free_if_valid ("free_sequence (twoBit.mBlockstarts)",     _seq->twoBit.mBlockstarts);
	free_if_valid ("free_sequence (twoBit.mBlocksizes)",      _seq->twoBit.mBlocksizes);

	if (_seq->hsx.fileInfo != NULL)
		{
		u32 fileNum;
		for (fileNum=0 ; fileNum<_seq->hsx.numFiles ; fileNum++)
			fclose_if_valid (_seq->hsx.fileInfo[fileNum].f);
		free_if_valid ("free_sequence (hsx.fileInfo)",        _seq->hsx.fileInfo);
		}

	                               free_if_valid ("free_sequence (partition.p)",    _seq->partition.p);
	if (_seq->partition.poolOwner) free_if_valid ("free_sequence (partition.pool)", _seq->partition.pool);

	free_if_valid ("free_sequence (_seq)", _seq);
	}

//----------
//
// load_sequence--
//	Load the next sequence from the associated file.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//
// Returns:
//	true if there was another sequence to load, false if not;  failures other
//	than eof result in fatality.
//
//----------

int load_sequence
   (seq*			_seq)
	{
	seqpartition*	sp;
	partition*		p;
	unspos			sepPos;
	char*			header;
	int				headerLen;
	unspos			oldTrueLen;

	debugNamesFile_1;

	if (_seq == NULL)                      suicide ("load_sequence(NULL)");
	if (_seq->preLoaded)                   { _seq->preLoaded = false;  return true; }
	if (_seq->fileType == seq_type_nofile) return false;
	if (!another_sequence (_seq))          return false;

	// get rid of sequence data from previous load

	_seq->len     = 0;
	_seq->trueLen = 0;

	//////////
	// read the sequence data, either
	//	o sequence is not partitioned => read just one sequence from the file
	//	o sequence is partitioned     => read all sequences from the file
	//////////

	sp = &_seq->partition;
	if ((sp->p == NULL) || (sp->state != seqpart_empty))
	 	{
		// not partitioned, just read the next sequence

		if (_seq->subsampleN > 0)
			{
			if (_seq->subsampleSkip > 0)
				skip_sequences (_seq, _seq->subsampleSkip);
			_seq->subsampleSkip = _seq->subsampleN-1;
			}

		load_sequence_core (_seq, /*keeper*/ true);
		}

	else
		{
		// partitioned, read the all sequences

		sp->state = seqpart_loading;

		// write first separator

		sequence_long_enough (_seq, 1, false);
		_seq->v[0] = 0;
		_seq->len  = 0;

		// read all the sequences

		while (another_sequence (_seq))
			{
			// load the next sequence

			debugPartitions_1;

			sepPos = _seq->len++;		// (advance length past the separator)
			oldTrueLen = _seq->trueLen;

			if (_seq->subsampleN > 0)
				{
				if (_seq->subsampleSkip > 0)
					skip_sequences (_seq, _seq->subsampleSkip);
				_seq->subsampleSkip = _seq->subsampleN-1;
				}
			load_sequence_core (_seq, /*keeper*/ true);

			header = (_seq->useFullNames)? _seq->header : _seq->shortHeader;
			headerLen = strlen(header);

			// add an entry to the partition table

			enough_partitions (_seq, sp->len+1, sp->poolLen+headerLen+1,
                               /*anticipate*/ true);

			p = &sp->p[sp->len];
			p->sepPos  = sepPos; 
			p->contig  = _seq->contig; 
			p->trueLen = _seq->trueLen - oldTrueLen;
			p->header  = sp->poolLen;
			strcpy (/*to*/ &sp->pool[p->header], /*from*/ header);
			sp->poolLen += headerLen+1;
			sp->len++;
			}

		// add final separator to table

		p = &sp->p[sp->len];
		p->sepPos = _seq->len; 

		sp->state = seqpart_ready;
		}

	// apply any required operators to it;  note that nib and 2bit sequences
	// are unmasked (if desired) during the earlier call to load_nib_sequence
	// or load_2bit_sequence, so there is no need to unmask them here;  further,
	// csfasta files do not have the concept of masking

	if ((_seq->doUnmask)
	 && ((_seq->fileType == seq_type_fasta)
	  || (_seq->fileType == seq_type_hsx)))
		upper_sequence (_seq);

	if (_seq->fileType == seq_type_qdna)
		{
		if ((_seq->softMaskFileName != NULL)
		 || (_seq->xMaskFileName    != NULL)
		 || (_seq->nMaskFileName    != NULL))
			suicidef ("masking not allowed for %s", sequence_filename(_seq));
		if (_seq->doUnmask)
			suicidef ("unmasking not allowed for %s", sequence_filename(_seq));
		if (((_seq->doRevCompFlags & rcf_comp) != 0)
		  && (_seq->qToComplement == NULL))
			suicidef ("reverse complement not allowed for %s\n",
			          "(the score file lacks complements)",
			           sequence_filename(_seq));
		}

	if (_seq->softMaskFileName != NULL)
		{
		if (_seq->softMaskComplement)
			mask_sequence_keep (_seq, _seq->softMaskFileName, -1);
		else
			mask_sequence      (_seq, _seq->softMaskFileName, -1);
		}
	if (_seq->xMaskFileName != NULL)
		{
		if (_seq->xMaskComplement)
			mask_sequence_keep (_seq, _seq->xMaskFileName, 'X');
		else
			mask_sequence      (_seq, _seq->xMaskFileName, 'X');
		}
	if (_seq->nMaskFileName != NULL)
		{
		if (_seq->nMaskComplement)
			mask_sequence_keep (_seq, _seq->nMaskFileName, 'N');
		else
			mask_sequence      (_seq, _seq->nMaskFileName, 'N');
		}

	if (_seq->doRevCompFlags == rcf_revcomp)
		rev_comp_sequence (_seq, _seq->qToComplement);
	else if (_seq->doRevCompFlags == rcf_rev)
		backward_sequence (_seq);
	else if (_seq->doRevCompFlags == rcf_comp)
		{
		backward_sequence (_seq);
		rev_comp_sequence (_seq, _seq->qToComplement);
		}

	debugPartitions_2;

	if (sequences_dbgDumpSequence)
		dump_sequence (stderr, _seq);

	return true;
	}


//-- skip_sequences--

static void skip_sequences
   (seq*	_seq,
	int		skipCount)
	{
	while ((skipCount-- > 0) && another_sequence_core (_seq))
		load_sequence_core (_seq, /*keeper*/ false);
	}


//-- load_sequence_core --

static void load_sequence_core
   (seq*	_seq,
	int		keeper)
	{
	// get rid of header data from previous load

	if (!_seq->lockedHeader)
		{
		if ((_seq->header != NULL) && (_seq->headerSize != 0))
			_seq->header[0] = 0;
		if ((_seq->shortHeader != NULL) && (_seq->shortHeaderSize != 0))
			_seq->shortHeader[0] = 0;
		}

	// read the next sequence for this type

	_seq->revCompFlags = rcf_forward;
	_seq->contig++;

	switch (_seq->fileType)
		{
		case seq_type_fasta:
			load_fasta_sequence (_seq, keeper);
			break;

		case seq_type_csfasta:
			load_csfasta_sequence (_seq, keeper);
			break;

		case seq_type_nib:
			load_nib_sequence (_seq, keeper);
			break;

		case seq_type_2bit:
			if (_seq->contigOfInterest != NULL)
				_seq->contig--; // (cancel earlier increment)
			load_2bit_sequence (_seq, keeper);
			break;

		case seq_type_hsx:
			if (_seq->contigOfInterest != NULL)
				_seq->contig--; // (cancel earlier increment)
			load_hsx_sequence (_seq, keeper);
			break;

		case seq_type_qdna:
			load_qdna_sequence (_seq, keeper);
			break;

		default:
			suicidef ("unknown sequence type: %X", _seq->fileType);
		}

	debugNamesFile_2;

	_seq->contigPending = false;

	if ((_seq->header != NULL)
	 && (_seq->headerSize != 0)
	 && ((_seq->shortHeader == NULL)
	  || (_seq->shortHeaderSize == 0)
	  || (_seq->shortHeader[0] == 0)
	  || (_seq->hasNickname)))
		create_short_header (_seq);
	}

//----------
//
// load_fasta_sequence--
//	Load the next fasta sequence from the associated file.
//
// A typical file looks like this:
//
//	> some header for the first sequence
//	GCGGTATCGCGCACAAGATTTAGGGATAGATCGTTTTGATGACCTCTCGCCACCTGGCAA
//	  ...
//	AAAAAAGGTAGGCCCATTAGCCCCCC
//
// The header line is optional.  However, if several sequences are included in
// the same file, the header lines are necessary as separators.  Nucleotides
// can be upper or lower case.  X can be used to indicate a masked position.
// whitespace can be added as desired (and so line breaks can be anywhere).
// digits are ignored so it is easy to use numerically annotated sequences.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//	int		keeper:	true  => actually load the sequence
//					false => just skip the sequence
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static void load_fasta_sequence
   (seq*	_seq,
	int		keeper)
	{
	u32		headerLen;
	int		ch;
	unspos	length;

	if (_seq == NULL) suicide ("load_fasta_sequence(NULL)");

	//////////
	// read the header
	//////////

	ch = skip_seq_whitespace (_seq);

	if (ch != '>')
		seq_ungetc (ch, _seq);
	else if (_seq->lockedHeader)
		{
		while ((ch != '\n') && (ch != '\r'))
			ch = seq_getc (_seq);
		}
	else
		{
		headerLen = 0;
		_seq->headerOwner = _seq->shortHeaderOwner = true;
		if (sequences_keepFastaArrow)
			{
			append_char (&_seq->header, &_seq->headerSize, &headerLen, '>');
			ch = seq_getc (_seq);
			while ((ch == ' ') || (ch == '\t'))
				{
				append_char (&_seq->header, &_seq->headerSize, &headerLen, ch);
				ch = seq_getc (_seq);
				}
			}
		else
			ch = skip_seq_whitespace (_seq);
		while ((ch != '\n') && (ch != '\r'))
			{
			append_char (&_seq->header, &_seq->headerSize, &headerLen, ch);
			ch = seq_getc (_seq);
			}
		append_char (&_seq->header, &_seq->headerSize, &headerLen, 0);

		if (_seq->nameTrigger != NULL)
			{
			char* triggerFound, *src, *dst;
			triggerFound = strstr (_seq->header, _seq->nameTrigger);
			if (triggerFound != NULL)
				{
				triggerFound += strlen (_seq->nameTrigger);
				for (src=triggerFound,dst=_seq->header ; *src!=0 ; )
					{
					ch = *(src++);
					if ((!isalnum(ch)) && (ch != '_'))
						break;
					*(dst++) = ch;
					}
				*dst = 0;
				headerLen = strlen (_seq->header);
				}
			}
		}

	//////////
	// read ahead to determine the length of the sequence and to pre-allocate
	// the vector
	//////////

	if (_seq->rewindable == -1)
		_seq->rewindable = (test_rewindability (_seq) == 0);

	if ((_seq->rewindable == true) && (keeper))
		{
		// read ahead, counting chars needed

		save_fstate (_seq);
		length = parse_fasta (_seq, /*storeEm*/ false);
		restore_fstate (_seq);

		// allocate the vector

		if ((length > maxSequenceLen) || (_seq->len > maxSequenceLen - length))
			suicidef ("in load_fasta_sequence for %s, "
			          "sequence length %s+%s exceeds maximum (%s)",
			          sequence_filename(_seq),
			          commatize(_seq->len), commatize(length),
			          commatize(maxSequenceLen));

		sequence_long_enough (_seq, _seq->len+length, false);
		}

	//////////
	// read the sequence
	//////////

	parse_fasta (_seq, /*storeEm*/ keeper);
	}

//----------
//
// parse_fasta--
//	Parse a fasta sequence from the associated file.  This assumes that the
//	sequence's file is positioned as the first character in the sequence,
//	*after* the sequence header line.
//
// (see load_fasta_sequence() for info about the file format)
//
//----------
//
// Arguments:
//	seq*	_seq:		The sequence to parse.
//	int		storeEm:	true  => store the results (in the sequence)
//						false => just count
//
// Returns:
//  The number of characters read into the sequence;  failure causes program
//	fatality.
//
//----------

static unspos parse_fasta
   (seq*	_seq,
	int		storeEm)
	{
	unspos	index, startLimit, endLimit, count;
	int		prevCh, ch;

	index      = 0;
	count      = 0;
	startLimit = _seq->startLimit;
	endLimit   = _seq->endLimit;

	// scan the file, keeping characters that are (a) nucleotides and (b) are
	// within our index limits

	prevCh = '\n';
	ch     = skip_fasta_whitespace (_seq);
	while (ch != EOF)
		{
		if ((prevCh == '\n') && (ch == '>')) // (start of next sequence)
			{ seq_ungetc (ch, _seq);  break; }

		switch (char_to_fasta_type[(u8)ch])
			{
			case _nucleotide:
				break;
			case _ambiguous:
				if (!_seq->allowAmbiDNA) goto bad_char;
				break;
			case _newline:
				ch = '\n'; // (allow for unix, mac, or pc line ends)
				goto next_char;
			case _bad:
				goto bad_char;
			}

		// this is a nucleotide, do we want it?

		index++;

		if ((startLimit != 0) && (index < startLimit)) goto next_char;
		if ((endLimit   != 0) && (index > endLimit))   goto next_char;

		// we want it;  are we just counting?

		if ((!storeEm) && (count+1 < count))
			suicidef ("in parse_fasta, "
			          "sequence length " unsposFmt "+1 overflows internal data type",
			          count);

		count++;
		if (!storeEm) goto next_char;

		// ok, let's store it

		if (_seq->len > maxSequenceLen - 1)
			suicidef ("in parse_fasta, "
			          "sequence length " unsposFmt "+1 exceeds maximum (" unsposFmt ")",
			          _seq->len, maxSequenceLen);

		sequence_long_enough (_seq, _seq->len+1, true);
		_seq->v[_seq->len++] = ch;

		// go try the next character

	next_char:
		prevCh = ch;
		ch = skip_fasta_whitespace (_seq);
		}

	if (storeEm)
		{
		_seq->v[_seq->len] = 0;				// (set the terminating zero)
		_seq->trueLen += index;				// (account for the characters
											//  .. we've read so far)
		}

	// make sure we got somethin' useful

	if ((startLimit != 0) && (startLimit > index))
		suicidef ("beyond end in %s (%08lX > %08lX)",
		          sequence_filename(_seq), startLimit, index);

	if ((endLimit != 0) && (endLimit > index))
		{
		if (_seq->endIsSoft)
			{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
		else
			suicidef ("beyond end in %s (%08lX > %08lX)",
			          sequence_filename(_seq), endLimit, index);
		}

	if ((count == 0) && (storeEm))
		{
		if (_seq->header == NULL)
			fprintf (stderr, "WARNING. %s contains an empty sequence\n",
			                 sequence_filename(_seq));
		else
			fprintf (stderr, "WARNING. %s contains an empty sequence:\n%s\n",
			                 sequence_filename(_seq), _seq->header);
		}

	if (startLimit == 0) _seq->start = 1;
	                else _seq->start = startLimit;

	// skip to the next sequence

	if ((storeEm) && (_seq->needTrueLen))
		{
		prevCh = '\n';
		ch     = skip_fasta_whitespace (_seq);
		while (ch != EOF)
			{
			if ((prevCh == '\n') && (ch == '>')) // (start of next sequence)
				{ seq_ungetc (ch, _seq);  break; }

			switch (char_to_fasta_type[(u8)ch])
				{
				case _nucleotide:
				case _ambiguous:
					_seq->trueLen++;
					break;
				case _newline:
					ch = '\n'; // (allow for unix, mac, or pc line ends)
					break;
				case _bad:
					goto bad_char;
				}

			// go try the next character

			prevCh = ch;
			ch = skip_fasta_whitespace (_seq);
			}
		}

	return count;

	// failure exits
	// $$$ report line number here

bad_char:
	if ((_seq->header == NULL) || (_seq->header[0] == 0))
		{
		if (dna_isprint(ch))
			suicidef ("bad fasta character in %s: %c",
				  sequence_filename(_seq), (int) ch);
		else
			suicidef ("bad fasta character in %s (ascii %02X)",
					  sequence_filename(_seq), (u8) ch);
		}
	else
		{
		if (dna_isprint(ch))
			suicidef ("bad fasta character in %s, %s: %c",
				  sequence_filename(_seq), _seq->header, (int) ch);
		else
			suicidef ("bad fasta character in %s, %s (ascii %02X)",
					  sequence_filename(_seq), _seq->header, (u8) ch);
		}

	return 0; // (never gets here)
	}

//----------
//
// load_csfasta_sequence--
//	Load the next fasta color sequence from the associated file.
//
// A typical file looks like this:
//
//	# Wed Apr 22 15:07:58 2009 ...
//	>538_743_229_F7
//	T013131021212033022020113200231003030002
//	>538_4021_559_F7
//	T002120310210323111000110101233231231210
//	>534_6488_139_F7
//	T112211320333111020130303120302210313113
//
// Line beginning with '#' are comments and are ignored, but they can only
// occur immediately in front of a header line. Lines beginning with ">" are
// header lines.  If the file contains only one sequence, the header line is
// optional.  Sequences must begin with a nucleotide and thereafter consist
// only of the digits '0', '1', '2' and '3'.  Sequences may occupy multiple
// lines.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//	int		keeper:	true  => actually load the sequence
//					false => just skip the sequence
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static void load_csfasta_sequence
   (seq*	_seq,
	int		keeper)
	{
	u32		headerLen;
	int		ch;
	unspos	length;

	if (_seq == NULL) suicide ("load_csfasta_sequence(NULL)");

	//////////
	// read the header
	//////////

	while (true) // skip comment lines
		{
		ch = skip_seq_whitespace (_seq);
		if (ch != '#') break;

		while ((ch != '\n') && (ch != '\r'))
			ch = seq_getc (_seq);
		}

	if (ch != '>')
		seq_ungetc (ch, _seq);
	else if (_seq->lockedHeader)
		{
		while ((ch != '\n') && (ch != '\r'))
			ch = seq_getc (_seq);
		}
	else
		{
		headerLen = 0;
		_seq->headerOwner = _seq->shortHeaderOwner = true;
		if (sequences_keepFastaArrow)
			{
			append_char (&_seq->header, &_seq->headerSize, &headerLen, '>');
			ch = seq_getc (_seq);
			while ((ch == ' ') || (ch == '\t'))
				{
				append_char (&_seq->header, &_seq->headerSize, &headerLen, ch);
				ch = seq_getc (_seq);
				}
			}
		else
			ch = skip_seq_whitespace (_seq);
		while ((ch != '\n') && (ch != '\r'))
			{
			append_char (&_seq->header, &_seq->headerSize, &headerLen, ch);
			ch = seq_getc (_seq);
			}
		append_char (&_seq->header, &_seq->headerSize, &headerLen, 0);

		if (_seq->nameTrigger != NULL)
			{
			char* triggerFound, *src, *dst;
			triggerFound = strstr (_seq->header, _seq->nameTrigger);
			if (triggerFound != NULL)
				{
				triggerFound += strlen (_seq->nameTrigger);
				for (src=triggerFound,dst=_seq->header ; *src!=0 ; )
					{
					ch = *(src++);
					if ((!isalnum(ch)) && (ch != '_'))
						break;
					*(dst++) = ch;
					}
				*dst = 0;
				headerLen = strlen (_seq->header);
				}
			}
		}

	//////////
	// read ahead to determine the length of the sequence and to pre-allocate
	// the vector
	//////////

	if (_seq->rewindable == -1)
		_seq->rewindable = (test_rewindability (_seq) == 0);

	if ((_seq->rewindable == true) && (keeper))
		{
		// read ahead, counting chars needed

		save_fstate (_seq);
		length = parse_csfasta (_seq, /*storeEm*/ false);
		restore_fstate (_seq);

		// allocate the vector

		if ((length > maxSequenceLen) || (_seq->len > maxSequenceLen - length))
			suicidef ("in load_csfasta_sequence for %s, "
			          "sequence length %s+%s exceeds maximum (%s)",
			          sequence_filename(_seq),
			          commatize(_seq->len), commatize(length),
			          commatize(maxSequenceLen));

		sequence_long_enough (_seq, _seq->len+length, false);
		}

	//////////
	// read the sequence
	//////////

	parse_csfasta (_seq, /*storeEm*/ keeper);
	}

//----------
//
// parse_csfasta--
//	Parse a csfasta sequence from the associated file.  This assumes that the
//	sequence's file is positioned as the first character in the sequence,
//	*after* the sequence header line.
//
// (see load_csfasta_sequence() for info about the file format)
//
//----------
//
// Arguments:
//	seq*	_seq:		The sequence to parse.
//	int		storeEm:	true  => store the results (in the sequence)
//						false => just count
//
// Returns:
//  The number of characters read into the sequence;  failure causes program
//	fatality.
//
//----------

static unspos parse_csfasta
   (seq*	_seq,
	int		storeEm)
	{
	unspos	index, startLimit, endLimit, count;
	int		prevCh, ch;
	u8		chType;

	index      = 0;
	count      = 0;
	startLimit = _seq->startLimit;
	endLimit   = _seq->endLimit;

	// scan the file, keeping characters that are (a) colors (or an initial
	// nucleotide) and (b) are within our index limits

	prevCh = '\n';
	ch     = skip_fasta_whitespace (_seq);
	while (ch != EOF)
		{
		if ((prevCh == '\n')
		 && ((ch == '#') || (ch == '>'))) 			// (start of next sequence)
			{ seq_ungetc (ch, _seq);  break; }

		chType = char_to_csfasta_type[(u8)ch];
		switch (chType)
			{
			case _nucleotide:
			case _color:
				break;
			case _newline:
				ch = '\n'; // (allow for unix, mac, or pc line ends)
				goto next_char;
			case _bad:
				goto bad_char;
			}

		// this is a color or nucleotide, do we want it?

		if ((index == 0) != (chType == _nucleotide))
			{
			if (index == 0) goto bad_nucleotide;
			           else goto bad_color;
			}

		index++;

		if ((startLimit != 0) && (index < startLimit)) goto next_char;
		if ((endLimit   != 0) && (index > endLimit))   goto next_char;

		// we want it;  are we just counting?

		if ((!storeEm) && (count+1 < count))
			suicidef ("in parse_csfasta, "
			          "sequence length " unsposFmt "+1 overflows internal data type",
			          count);

		count++;
		if (!storeEm) goto next_char;

		// ok, let's store it

		if (_seq->len > maxSequenceLen - 1)
			suicidef ("in parse_csfasta, "
			          "sequence length " unsposFmt "+1 exceeds maximum (" unsposFmt ")",
			          _seq->len, maxSequenceLen);

		sequence_long_enough (_seq, _seq->len+1, true);
		_seq->v[_seq->len++] = ch;

		// go try the next character

	next_char:
		prevCh = ch;
		ch = skip_fasta_whitespace (_seq);
		}

	if (storeEm)
		{
		_seq->v[_seq->len] = 0;				// (set the terminating zero)
		_seq->trueLen += index;				// (account for the characters
											//  .. we've read so far)
		}

	// make sure we got somethin' useful

	if ((startLimit != 0) && (startLimit > index))
		suicidef ("beyond end in %s (%08lX > %08lX)",
		          sequence_filename(_seq), startLimit, index);

	if ((endLimit != 0) && (endLimit > index))
		{
		if (_seq->endIsSoft)
			{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
		else
			suicidef ("beyond end in %s (%08lX > %08lX)",
			          sequence_filename(_seq), endLimit, index);
		}

	if ((count == 0) && (storeEm))
		{
		if (_seq->header == NULL)
			fprintf (stderr, "WARNING. %s contains an empty sequence\n",
			                 sequence_filename(_seq));
		else
			fprintf (stderr, "WARNING. %s contains an empty sequence:\n%s\n",
			                 sequence_filename(_seq), _seq->header);
		}

	if (startLimit == 0) _seq->start = 1;
	                else _seq->start = startLimit;

	// skip to the next sequence

	if ((storeEm) && (_seq->needTrueLen))
		{
		prevCh = '\n';
		ch     = skip_fasta_whitespace (_seq);
		while (ch != EOF)
			{
			if ((prevCh == '\n')
			 && ((ch == '#') || (ch == '>')))		// (start of next sequence)
				{ seq_ungetc (ch, _seq);  break; }

			switch (char_to_csfasta_type[(u8)ch])
				{
				case _nucleotide:
					goto bad_color;
				case _color:
					_seq->trueLen++;
					break;
				case _newline:
					ch = '\n'; // (allow for unix, mac, or pc line ends)
					break;
				case _bad:
					goto bad_char;
				}

			// go try the next character

			prevCh = ch;
			ch = skip_fasta_whitespace (_seq);
			}
		}

	return count;

	// failure exits
	// $$$ report line number and sequence name here

bad_char:
	if (dna_isprint(ch))
		suicidef ("bad csfasta character in %s: %c",
		          sequence_filename(_seq), (int) ch);
	else
		suicidef ("bad csfasta character in %s (ascii %02X)",
		          sequence_filename(_seq), (u8) ch);
	return 0; // (never gets here)

bad_nucleotide:
	if (dna_isprint(ch))
		suicidef ("bad csfasta nucleotide in %s: %c",
		          sequence_filename(_seq), (u8) ch);
	else
		suicidef ("bad csfasta nucleotide in %s (ascii %02X)",
		          sequence_filename(_seq), (int) ch);
	return 0; // (never gets here)

bad_color:
	if (dna_isprint(ch))
		suicidef ("bad csfasta color in %s: %c",
		          sequence_filename(_seq), (int) ch);
	else
		suicidef ("bad csfasta color in %s (ascii %02X)",
		          sequence_filename(_seq), (u8) ch);
	return 0; // (never gets here)
	}

//----------
//
// load_nib_sequence--
//	Load a nib sequence from the associated file.
//
// A nib file stores each nucleotide in four bits (one nybble).  The file
// consists of a 4 byte magic number, followed by a 4 byte length, followed by
// the nucleotides.  The magic number is in the file as
//	(first byte) 3A 3D E9 6B (third byte)
// The length field is in little-endian order, so
//	(first byte) C0 E1 E4 00 (third byte)
// means 0x00E4E1C0 bytes (15 million).  The length is the number of
// nucleotides.  The first nucleotide is in the most significant nybble of the
// 9th byte, the second one is in the least significant nybble, the third in
// the 10th byte (msnybble), and so on.  Nybble bits are mapped to characters
// as per the tables nibTo1stChar[] and nibTo2ndChar[].
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//	int		keeper:	true  => actually load the sequence
//					false => just skip the sequence
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static void load_nib_sequence
   (seq*		_seq,
	int			keeper)
	{
	u32			magic, length;
	unspos		newSeqLen, ix;
	u32			startLimit, endLimit, startIndex;
	u32			bytesLeft, bytesToRead, bytesRead;
	u8			ch;
	const u8*	to1stChar, *to2ndChar;
	u8*			dst;

	if (!keeper) return;
	if (_seq == NULL) suicide ("load_nib_sequence(NULL)");

	//////////
	// get the sequence length
	//////////

	// check the magic number and decide if it's little or big endian

	magic = read_4_big (_seq);

	// read the length

	if (magic == nibMagicLittle)
		length = read_4_little (_seq);
	else if (magic == nibMagicBig)
		length = read_4_big (_seq);
	else
		{
		length = 0; // (placate compiler)
		suicidef ("bad nib magic number in %s (%08lX)",
		          sequence_filename(_seq), magic);
		}

	if ((length == 0) || (((s32) length) == -1))
		suicidef ("bad nib length in %s (%08lX)", sequence_filename(_seq), length);

	// validate sequence limits

	if ((_seq->startLimit != 0) && (_seq->startLimit > (unspos) length))
		suicidef ("beyond end in %s (%ld > %ld)",
		          sequence_filename(_seq), _seq->startLimit, length);

	if ((_seq->endLimit != 0) && (_seq->endLimit > (unspos) length))
		{
		if (_seq->endIsSoft)
			{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
		else
			suicidef ("beyond end in %s (%ld > %ld)",
			          sequence_filename(_seq), _seq->endLimit, length);
		}

	startLimit = (u32) _seq->startLimit;
	endLimit   = (u32) _seq->endLimit;

	_seq->trueLen += length;

	// skip ahead to the first desired base, and determine how many bases
	// we'll read

	if (startLimit == 0) startLimit = 1;
	startIndex = startLimit - 1;

	bytesLeft = length;
	if (startIndex > 0)
		{ skip_chars (_seq, startIndex/2);  bytesLeft -= 2*(startIndex/2); }

	length = bytesLeft;
	if ((startIndex&1) != 0)	// start offset is odd
		length--;
	if (endLimit != 0)
		{
		if (length > endLimit - startIndex)
			length = endLimit - startIndex;
		}

	//////////
	// allocate the vector, including an extra byte since we may overshoot by
	// 1 when unpacking
	//////////

#if (maxSequenceIndex <= 32)	// otherwise compiler complains that this test is
								// .. always false
	if ((length > maxSequenceLen) || (_seq->len > maxSequenceLen - length))
		suicidef ("in load_nib_sequence for %s, "
		          "sequence length %s+%s exceeds maximum (%s)",
		          sequence_filename(_seq),
		          commatize(_seq->len), commatize(length),
		          commatize(maxSequenceLen));
#endif

	newSeqLen = _seq->len + length;
	sequence_long_enough (_seq, newSeqLen+1, false);

	//////////
	// read the sequence
	//////////

	// decide which lookup tables we'll use

	if (_seq->doUnmask)
		{
		to1stChar = nibTo1stCharUnmasked;
		to2ndChar = nibTo2ndCharUnmasked;
		}
	else
		{
		to1stChar = nibTo1stChar;
		to2ndChar = nibTo2ndChar;
		}

	// read the first, partial, byte

	ix = _seq->len;
	if ((startIndex&1) != 0)	// start offset is odd
		{
		ch = seq_getc (_seq);  bytesLeft -= 2;
		_seq->v[ix++] = to2ndChar[ch];
		}

	// process any bytes in the pending buffer, one at a time

	while ((ix < newSeqLen) && (_seq->pendingLen > 0))
		{
		ch = seq_getc (_seq);  bytesLeft -= 2;
		_seq->v[ix++] = to1stChar[ch];
		_seq->v[ix++] = to2ndChar[ch];
		}

	// read the remaining bytes to the tail end of the buffer

	bytesToRead = ((newSeqLen-ix) + 1) / 2;
	dst = _seq->v + _seq->size - bytesToRead;
	if (bytesToRead > 0)
		{
		bytesRead = fread (dst, 1, bytesToRead, _seq->f);
		if (bytesRead != bytesToRead)
			suicidef ("in load_nib_sequence(%s), block read\n"
			          "wanted %d bytes, only got %d",
			          sequence_filename(_seq), bytesToRead, bytesRead);
		}

	// unpack those bytes;  note that although we are writing into the same
	// buffer that we are reading from, and writing two bytes for each one read,
	// the write pointer will not overtake the read pointer;  further note that
	// we may unpack an extra nybble (this will be overwritten when we set the
	// terminating zero)

	while (ix < newSeqLen)
		{
		ch = *(dst++);  bytesLeft--;
		_seq->v[ix++] = to1stChar[ch];
		_seq->v[ix++] = to2ndChar[ch];
		}

	_seq->v[newSeqLen] = 0;					// (set the terminating zero)
	_seq->len = newSeqLen;

	skip_chars (_seq, bytesLeft);			// skip the rest of the sequence
	_seq->pendingLen   = 0;					// (discard any pending chars)
	_seq->pendingStack = _seq->pendingChars + seqBufferSize;

	if (startLimit == 0) _seq->start = 1;
	                else _seq->start = startLimit;

	//////////
	// create a header
	//////////

	if (!_seq->lockedHeader)
		{
		length = snprintf (_seq->header, 0, "%s:" unsposDashFmt,
		                   sequence_filename(_seq), _seq->start, _seq->start + _seq->len-1);

		if (_seq->headerSize < length+1)
			{
			_seq->header      = realloc_or_die ("load_nib_sequence (header)",
			                                   _seq->header, length+1);
			_seq->headerSize  = length+1;
			}
		_seq->headerOwner = _seq->shortHeaderOwner = true;

		snprintf (_seq->header, length+1, "%s:" unsposDashFmt,
		          sequence_filename(_seq), _seq->start, _seq->start + _seq->len-1);
		}
	}

//----------
//
// read_2bit_header, load_2bit_sequence--
//	Load a 2bit sequence from the associated file.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//	int		keeper:	(load_2bit_sequence only)
//					true  => actually load the sequence
//					false => just skip the sequence
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static int find_2bit_sequence    (seq* _seq, char* name);
static u32 read_2bit_index_entry (seq* _seq, char seqName[256], u32 seqNum);

//--- read_2bit_header ---

static void read_2bit_header
   (seq*	_seq)
	{
	u32		magic, version, reserved;

	// read and validate the header

	magic = read_4_big (_seq);

	if (magic == twobitMagicLittle)
		_seq->twoBit.bigEndian = false;
	else if (magic == twobitMagicBig)
		_seq->twoBit.bigEndian = true;
	else
		suicidef ("bad 2bit magic number in %s (%08lX)",
				  sequence_filename(_seq), magic);

	version                 = read_4 (_seq, _seq->twoBit.bigEndian);
	_seq->twoBit.numContigs = read_4 (_seq, _seq->twoBit.bigEndian);
	reserved                = read_4 (_seq, _seq->twoBit.bigEndian);

	if (version != 0)
		suicidef ("bad 2bit version in %s (%08lX)",
		          sequence_filename(_seq), version);
	if (reserved != 0)
		suicidef ("bad 2bit header word 4 in %s (%08lX)",
		          sequence_filename(_seq), reserved);

	if (_seq->twoBit.numContigs == 0)
		suicidef ("empty 2bit file %s", sequence_filename(_seq));

	// save index's file position

	_seq->twoBit.indexFilePos = _seq->twoBit.contigFilePos = ftell (_seq->f);

	// if we have a single contig-of-interest, locate it

	if (_seq->contigOfInterest != NULL)
		{
		if (!find_2bit_sequence (_seq, _seq->contigOfInterest))
			suicidef ("2bit file %s doesn't contain %s",
					  sequence_filename(_seq), _seq->contigOfInterest);
		}

	}

//--- load_2bit_sequence ---

static void load_2bit_sequence
   (seq*	_seq,
	int		keeper)
	{
	char	seqName[maxSequenceName+1];
	int		numChars;
	u32		dnaSize, reserved;
	u32		nBlockCount, maskBlockCount;
	u32		seqDataPos, seekPos;
	unspos	oldSeqLen, ix;
	u32		startLimit, endLimit, length, startIndex, endIndex;
	u32		basesToGo, bytesToSkip, bytesToRead, bytesRead;
	u32		blockIx, s, e, scanIx;
	u8		ch;
	u8*		data, *dst;
	char*	seekType;
	int		err;

	_seq->pendingLen   = 0;					// (discard any pending chars)
	_seq->pendingStack = _seq->pendingChars + seqBufferSize;

	//////////
	// read the sequence's index table entry
	//////////

	err = fseek (_seq->f, _seq->twoBit.contigFilePos, SEEK_SET);
	if (err != 0)
		{ seekType = "index";  seekPos = _seq->twoBit.contigFilePos;  goto fseek_failed; }

	seqDataPos = read_2bit_index_entry (_seq, seqName, _seq->contig);
	_seq->twoBit.contigFilePos = ftell (_seq->f);

	if (!keeper) return;

	// copy the sequence name as our header

	numChars = strlen (seqName);
	if (_seq->headerSize < (unsigned) (numChars+1))
		{
		_seq->header      = realloc_or_die ("load_2bit_sequence (header)",
						 				  _seq->header, numChars+1);
		_seq->headerSize  = numChars+1;
		_seq->headerOwner = _seq->shortHeaderOwner = true;
		}

	strcpy (/*to*/ _seq->header, /*from*/ seqName);

	//////////
	// make sure we have enough room for the sequence's data
	//////////

	err = fseek (_seq->f, seqDataPos, SEEK_SET);
	if (err != 0)
		{ seekType = "header data";  seekPos = seqDataPos;  goto fseek_failed; }

	dnaSize = read_4 (_seq, _seq->twoBit.bigEndian);
	seqDataPos += 4;

	if ((dnaSize == 0) || (((s32) dnaSize) == -1))
		suicidef ("bad 2bit length in %s (%08lX)", sequence_filename(_seq), dnaSize);

	if ((_seq->startLimit != 0) && (_seq->startLimit > (unspos) dnaSize))
		suicidef ("beyond end in %s (%ld > %ld)",
		          sequence_filename(_seq), _seq->startLimit, dnaSize);

	if ((_seq->endLimit != 0) && (_seq->endLimit > (unspos) dnaSize))
		{
		if (_seq->endIsSoft)
			{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
		else
			suicidef ("beyond end in %s (%ld > %ld)",
			          sequence_filename(_seq), _seq->endLimit, dnaSize);
		}

	startLimit = (u32) _seq->startLimit;
	endLimit   = (u32) _seq->endLimit;

	if (startLimit == 0) startLimit = 1;
	if (endLimit   == 0) endLimit   = dnaSize;

	_seq->trueLen += dnaSize;

	// allocate the vector; we ask for an additional three characters because
	// during unpacking we write 4 bytes at a time, and thus may overshoot the
	// end by 3
	// $$$ note that this may be way more than needed, if start and end limits
	// $$$ .. reduce the number of characters we actually want;  this could be
	// $$$ .. improved if it causes problems

#if (maxSequenceIndex <= 32)	// otherwise compiler complains that this test is
								// .. always false
	if ((dnaSize+3 > maxSequenceLen) || (_seq->len > maxSequenceLen - (dnaSize+3)))
		goto sequence_too_big;
#endif

	sequence_long_enough (_seq, _seq->len + dnaSize+3, false);

	//////////
	// read and save the intervening block-marking fields
	//////////

	// make sure we have enough room for the n-blocks

	nBlockCount = read_4 (_seq, _seq->twoBit.bigEndian);
	seqDataPos += 4;

	if (nBlockCount > _seq->twoBit.nBlocksSize)
		{
		_seq->twoBit.nBlockStarts = (u32*) realloc_or_die ("nBlockStarts", _seq->twoBit.nBlockStarts, nBlockCount * sizeof(u32));
		_seq->twoBit.nBlockSizes  = (u32*) realloc_or_die ("nBlockSizes",  _seq->twoBit.nBlockSizes,  nBlockCount * sizeof(u32));
		_seq->twoBit.nBlocksSize  = nBlockCount;
		}

	// read the n-blocks

	for (blockIx=0 ; blockIx<nBlockCount ; blockIx++)
		_seq->twoBit.nBlockStarts[blockIx] = read_4 (_seq, _seq->twoBit.bigEndian);

	for (blockIx=0 ; blockIx<nBlockCount ; blockIx++)
		_seq->twoBit.nBlockSizes[blockIx] = read_4 (_seq, _seq->twoBit.bigEndian);

	seqDataPos += 4 * (2 * nBlockCount);

	// make sure we have enough room for the mask-blocks;  note that if we are
	// unmasking the sequence, then we skip over the mask-blocks

	maskBlockCount = read_4 (_seq, _seq->twoBit.bigEndian);
	seqDataPos += 4;

	if (_seq->doUnmask)
		{
		seqDataPos += 4 * (2 * maskBlockCount);
		err = fseek (_seq->f, seqDataPos, SEEK_SET);
		if (err != 0)
			{ seekType = "mask data";  seekPos = seqDataPos;  goto fseek_failed; }
		maskBlockCount = 0;
		}
	else
		{
		if (maskBlockCount > _seq->twoBit.mBlocksSize)
			{
			_seq->twoBit.mBlockstarts = (u32*) realloc_or_die ("mBlockstarts", _seq->twoBit.mBlockstarts, maskBlockCount * sizeof(u32));
			_seq->twoBit.mBlocksizes  = (u32*) realloc_or_die ("mBlocksizes",  _seq->twoBit.mBlocksizes,  maskBlockCount * sizeof(u32));
			_seq->twoBit.mBlocksSize  = maskBlockCount;
			}
		}

	// read the mask-blocks

	for (blockIx=0 ; blockIx<maskBlockCount ; blockIx++)
		_seq->twoBit.mBlockstarts[blockIx] = read_4 (_seq, _seq->twoBit.bigEndian);

	for (blockIx=0 ; blockIx<maskBlockCount ; blockIx++)
		_seq->twoBit.mBlocksizes[blockIx] = read_4 (_seq, _seq->twoBit.bigEndian);

	seqDataPos += 4 * (2 * maskBlockCount);

	// skip the reserved data prefix

	reserved = read_4 (_seq, _seq->twoBit.bigEndian);
	if (reserved != 0)
		suicidef ("bad 2bit reserved data prefix in %s\n"
		          "         (data at %08lX is %08lX)",
		          sequence_filename(_seq), seqDataPos, reserved);

	//////////
	// read the sequence's data
	//////////

	// skip to the first byte containing data of interest

	startIndex = startLimit-1;
	length     = basesToGo = endLimit+1 - startLimit;

	bytesToSkip = startIndex / 4;
	if (bytesToSkip != 0)
		{
		err = fseek (_seq->f, bytesToSkip, SEEK_CUR);
		if (err != 0)
			{ seekType = "data";  seekPos = _seq->twoBit.contigFilePos;  goto fseek_failed; }
		startIndex -= 4*bytesToSkip;
		}

	// read the leading partial byte (if any)

	ix = oldSeqLen = _seq->len;
	if (startIndex > 0)
		{
		ch = (u8) seq_getc (_seq);
		data = (u8*) (twobitToChars[ch] + startIndex);
		while (*data != 0)
			{
			if (basesToGo-- <= 0) break;
			_seq->v[ix++] = *(data++);
			}
		}

	// process any bytes in the pending buffer;  note that we may end up writing
	// as many as 3 bytes beyond the end of the sequence, but will correct this
	// when we write the terminating zero;  also note that we have to separate
	// the last iteration of this loop, since basesToGo is unsigned

	for ( ; (basesToGo>=4)&&(_seq->pendingLen>0) ; basesToGo-=4)
		{
		ch   = (u8)  seq_getc (_seq);
		data = (u8*) twobitToChars[ch];
		_seq->v[ix++] = data[0];
		_seq->v[ix++] = data[1];
		_seq->v[ix++] = data[2];
		_seq->v[ix++] = data[3];
		}

	if ((basesToGo > 0) && (_seq->pendingLen > 0))
		{
		ch   = (u8)  seq_getc (_seq);
		data = (u8*) twobitToChars[ch];
		_seq->v[ix++] = data[0];
		_seq->v[ix++] = data[1];
		_seq->v[ix++] = data[2];
		_seq->v[ix++] = data[3];
		basesToGo = 0;
		}

	// read the remaining bytes to the tail end of the buffer

	bytesToRead = (basesToGo + 3) / 4;
	dst = _seq->v + _seq->size - bytesToRead;
	if (bytesToRead > 0)
		{
		bytesRead = fread (dst, 1, bytesToRead, _seq->f);
		if (bytesRead != bytesToRead) goto read_failed;
		}

	// unpack those bytes;  note that although we are writing into the same
	// buffer that we are reading from, and writing four bytes for each one
	// read, the write pointer will not overtake the read pointer;  as above,
	// we may end up writing as many as 3 bytes beyond the end of the sequence

	for ( ; basesToGo>=4 ; basesToGo-=4)
		{
		ch   = *(dst++);
		data = (u8*) twobitToChars[ch];
		_seq->v[ix++] = data[0];
		_seq->v[ix++] = data[1];
		_seq->v[ix++] = data[2];
		_seq->v[ix++] = data[3];
		}

	if (basesToGo > 0)
		{
		ch   = *(dst++);
		data = (u8*) twobitToChars[ch];
		_seq->v[ix++] = data[0];
		_seq->v[ix++] = data[1];
		_seq->v[ix++] = data[2];
		_seq->v[ix++] = data[3];
		basesToGo = 0;
		}

	_seq->len += length;
	_seq->v[_seq->len] = 0;					// (set the terminating zero)

	//////////
	// mark the Ns and masked bases
	//////////

	startIndex = startLimit-1;
	endIndex   = endLimit;

	for (blockIx=0 ; blockIx<nBlockCount ; blockIx++)
		{
		s = _seq->twoBit.nBlockStarts[blockIx];
		e = s + _seq->twoBit.nBlockSizes[blockIx];
		if (e <= startIndex) continue;
		if (s >= endIndex)   continue;
		if (s <  startIndex) s = startIndex;
		if (e >  endIndex)   e = endIndex;
		s -= startIndex;
		e -= startIndex;
		for (scanIx=s ; scanIx<e ; scanIx++)
			_seq->v[oldSeqLen+scanIx] = 'N';
		}

	for (blockIx=0 ; blockIx<maskBlockCount ; blockIx++)
		{
		s = _seq->twoBit.mBlockstarts[blockIx];
		e = s + _seq->twoBit.mBlocksizes[blockIx];
		if (e <= startIndex) continue;
		if (s >= endIndex)   continue;
		if (s <  startIndex) s = startIndex;
		if (e >  endIndex)   e = endIndex;
		s -= startIndex;
		e -= startIndex;
		for (scanIx=s ; scanIx<e ; scanIx++)
			_seq->v[oldSeqLen+scanIx] = dna_tolower (_seq->v[oldSeqLen+scanIx]);
		}

	_seq->twoBit.contigLoaded = true;

	if (startLimit == 0) _seq->start = 1;
	                else _seq->start = startLimit;
	return;

// failure exits

fseek_failed:
	suicidef ("failed to seek to position in \"%s\"\n"
	          "in load_2bit_sequence, %s fseek(%08lX) returned %d",
	          sequence_filename(_seq), seekType, seekPos, err);
	return; // (never gets here)

sequence_too_big:
	suicidef ("in load_2bit_sequence for %s, "
			  "sequence length %s+%s exceeds maximum (%s)",
			  sequence_filename(_seq),
			  commatize(_seq->len),commatize(dnaSize+3),
			  commatize(maxSequenceLen));
	return; // (never gets here)

read_failed:
	suicidef ("in load_2bit_sequence for %s,"
	          " block read for sequence %u\n"
	          "wanted %d bytes, only got %d",
	          sequence_filename(_seq), _seq->contig, bytesToRead, bytesRead);
	return; // (never gets here)
	}

//--- find_2bit_sequence ---

static int find_2bit_sequence
   (seq*	_seq,
	char*	name)
	{
	char	seqName[maxSequenceName+1];
	u32		ix;

	for (ix=0 ; ix<_seq->twoBit.numContigs ; ix++)
		{
		_seq->twoBit.contigFilePos = ftell (_seq->f);
		read_2bit_index_entry (_seq, seqName, ix+1);
		if (strcmp (seqName, name) == 0)
			{ _seq->contig = ix + 1;  return true; }
		}

	return false;
	}

//--- read_2bit_index_entry ---

static u32 read_2bit_index_entry
   (seq*			_seq,
	char			seqName[maxSequenceName+1],
	u32				seqNum)
	{
	unsigned int	nameSize;
	size_t			bytesRead;

	// read the name

	nameSize = getc_or_die (_seq->f, _seq->fileName);
	if (nameSize > 0)
		{
		bytesRead = fread (seqName, 1, nameSize, _seq->f);
		if (bytesRead != nameSize)
			suicidef ("in load_2bit_sequence for %s, short read for sequence %u\n"
			          "wanted %d bytes, only got %d",
			          sequence_filename(_seq), seqNum, nameSize, bytesRead);
		}
	seqName[nameSize] = 0;

	// read the data offset

	return read_4 (_seq, _seq->twoBit.bigEndian);
	}

//----------
//
// read_hsx_header, load_hsx_sequence--
//	Load a sequence from the associated hsx file.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//	int		keeper:	(load_hsx_sequence only)
//					true  => actually load the sequence
//					false => just skip the sequence
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static u64   lookup_hsx_sequence  (seq* _seq, char* name);
static u64   find_hsx_sequence    (seq* _seq, char* name,
                                   u64 bucketStart, u64 bucketEnd);
static char* read_hsx_index_entry (seq* _seq);
static char* read_hsx_string      (seq* _seq, FILE* f);

//--- read_hsx_header ---

static void read_hsx_header
   (seq*		_seq)
	{
	u32			fileInfoOffset[255];
	u32			magic, headerLength;
	char*		s;
	char		extension[10];
	u64			fileOffset;
	long int	bucketOffset;
	u32			infoBytes, nameBytes;
	char*		slash, *dot, *nameScan;
	int  		baseLen, pathLen;
	u32			fileNum;
	int			err;

	// read and validate the header

	magic = read_4_big (_seq);

	if (magic == hsxMagicLittle)
		_seq->hsx.bigEndian = false;
	else if (magic == hsxMagicBig)
		_seq->hsx.bigEndian = true;
	else
		suicidef ("bad hsx magic number in %s (%08lX)",
				  sequence_filename(_seq), magic);

	_seq->hsx.version = read_4 (_seq, _seq->hsx.bigEndian);
	if (_seq->hsx.version != 0x00000100L)
		suicidef ("bad hsx version in %s (%08lX)",
		          sequence_filename(_seq), _seq->hsx.version);

	headerLength = read_4 (_seq, _seq->hsx.bigEndian);
	if (headerLength != 0x1C)
		suicidef ("bad hsx header length in %s (%08lX)",
		          sequence_filename(_seq), headerLength);

	_seq->hsx.numFiles        =       read_4 (_seq, _seq->hsx.bigEndian);
	_seq->hsx.fileTableOffset = (u64) read_4 (_seq, _seq->hsx.bigEndian);
	_seq->hsx.numBuckets      =       read_4 (_seq, _seq->hsx.bigEndian);
	_seq->hsx.hashTableOffset = (u64) read_4 (_seq, _seq->hsx.bigEndian);
	_seq->hsx.numContigs      =       read_4 (_seq, _seq->hsx.bigEndian);
	_seq->hsx.seqTableOffset  = (u64) read_4 (_seq, _seq->hsx.bigEndian);

	if (_seq->hsx.numFiles == 0)
		suicidef ("empty file table in hsx file %s", sequence_filename(_seq));

	if (_seq->hsx.numFiles > 255)
		suicidef ("corrupt header in hsx file %s (numFiles > 255;  %d)",
		          sequence_filename(_seq), _seq->hsx.numFiles);

	if (_seq->hsx.numBuckets == 0)
		suicidef ("corrupt header in hsx file %s (numBuckets = 0)",
		          sequence_filename(_seq));

	// read and validate the file table

	err = fseek (_seq->f, (long int) _seq->hsx.fileTableOffset, SEEK_SET);
	if (err != 0)
		suicidef ("in read_hsx_header for %s, file table fseek(%08lX) returned %d",
		          sequence_filename(_seq), _seq->hsx.fileTableOffset, err);

	for (fileNum=0 ; fileNum<_seq->hsx.numFiles ; fileNum++)
		fileInfoOffset[fileNum] = (u64) read_4 (_seq, _seq->hsx.bigEndian);

	slash = strrchr (_seq->fileName, pathSlash);
	dot   = strrchr (_seq->fileName, '.');
	if ((dot == NULL) || ((slash != NULL) && (dot < slash)))
		baseLen = strlen(_seq->fileName);
	else
		baseLen = dot - _seq->fileName;
	if (slash == NULL)
		pathLen = 0;
	else
		pathLen = slash+1 - _seq->fileName;

	infoBytes = sizeof(hsxfileinfo) * _seq->hsx.numFiles;
	nameBytes = 0;
	for (fileNum=0 ; fileNum<_seq->hsx.numFiles ; fileNum++)
		{
		err = fseek (_seq->f, (long int) fileInfoOffset[fileNum], SEEK_SET);
		if (err != 0)
			suicidef ("in read_hsx_header for %s, file table fseek(%08lX) returned %d",
			          sequence_filename(_seq), fileInfoOffset[fileNum], err);

		s = read_hsx_string (_seq, _seq->f);
		if ((strcmp (s, "fa")    != 0)
		 && (strcmp (s, "fasta") != 0))
			suicidef ("in read_hsx_header for %s, unsupported file type: %s",
					  sequence_filename(_seq), s);
		strncpy (/*to*/ extension, /*from*/ s, sizeof(extension));

		s = read_hsx_string (_seq, _seq->f);
		if (s[0] != 0)
			nameBytes += pathLen + strlen(s) + 1 + strlen(extension) + 1;
		else
			nameBytes += baseLen + 1 + strlen(extension) + 1;
		}

	_seq->hsx.fileInfo = (hsxfileinfo*) zalloc_or_die ("read_hsx_header", infoBytes + nameBytes);

	nameScan = ((char*) _seq->hsx.fileInfo) + infoBytes;
	for (fileNum=0 ; fileNum<_seq->hsx.numFiles ; fileNum++)
		{
		_seq->hsx.fileInfo[fileNum].name = nameScan;
		_seq->hsx.fileInfo[fileNum].f    = NULL;

		err = fseek (_seq->f, (long int) fileInfoOffset[fileNum], SEEK_SET);
		if (err != 0)
			suicidef ("in read_hsx_header for %s, file table fseek(%08lX) returned %d",
			          sequence_filename(_seq), fileInfoOffset[fileNum], err);

		s = read_hsx_string (_seq, _seq->f);
		strncpy (/*to*/ extension, /*from*/ s, sizeof(extension));

		s = read_hsx_string (_seq, _seq->f);
		if (s[0] != 0)
			{
			strncpy (/*to*/    nameScan,
					 /*from*/  _seq->fileName,
					 /*limit*/ pathLen);
			strcpy  (/*to*/   nameScan + pathLen,
					 /*from*/ s);
			nameScan[pathLen+strlen(s)] = '.';
			strcpy  (/*to*/   nameScan + pathLen+strlen(s) + 1,
					 /*from*/ extension);
			nameScan += pathLen + strlen(s) + 1 + strlen(extension) + 1;
			}
		else
			{
			strncpy (/*to*/    nameScan,
					 /*from*/  _seq->fileName,
					 /*limit*/ baseLen);
			nameScan[baseLen] = '.';
			strcpy  (/*to*/   nameScan + baseLen + 1,
					 /*from*/ extension);
			nameScan += baseLen + 1 + strlen(extension) + 1;
			}
		}

	// if we have a single contig-of-interest, locate it

	if (_seq->contigOfInterest != NULL)
		{
		fileOffset = lookup_hsx_sequence (_seq, _seq->contigOfInterest);
		if ((fileOffset & hsxMsBit5) != 0)
			suicidef ("hsx file %s doesn't contain %s",
					  sequence_filename(_seq), _seq->contigOfInterest);
		if (fileOffset > hsxMaxFilePos)
			suicidef ("in read_hsx_header for %s,"
			          " file pos for %s (%010lX) exceeds max (%010lX)",
					  sequence_filename(_seq), _seq->contigOfInterest,
					  fileOffset, hsxMaxFilePos);
		_seq->hsx.contigFilePos = fileOffset;
		}

	// otherwise, if we have no contig names, locate the first sequence in
	// the index

	else if (_seq->namesFileName == NULL)
		{
		bucketOffset = (long int) _seq->hsx.hashTableOffset;
		err = fseek (_seq->f, bucketOffset, SEEK_SET);
		if (err != 0)
			suicidef ("in read_hsx_header for %s,"
			          " file table fseek(%010lX) returned %d",
					  sequence_filename(_seq), 0, err);

		fileOffset = read_5 (_seq, _seq->hsx.bigEndian) & ~hsxMsBit5;
		if (fileOffset > hsxMaxFilePos)
			suicidef ("in read_hsx_header for %s,"
			          " file pos for index 0 (%010lX) exceeds max (%010lX)",
					  sequence_filename(_seq), fileOffset, hsxMaxFilePos);
		_seq->hsx.contigFilePos = fileOffset;
		}

	}

//--- load_hsx_sequence ---

static void load_hsx_sequence
   (seq*	_seq,
	int		keeper)
	{
	int		err;
	char*	seqName;
	int		numChars;
	unspos	dnaSize;
	unspos	index, startLimit, endLimit;
	int		prevCh, ch;
	char*	seqFName;
	FILE*	seqF;

	_seq->pendingLen   = 0;					// (discard any pending chars)
	_seq->pendingStack = _seq->pendingChars + seqBufferSize;

	//////////
	// read the sequence's index table entry
	//////////

	if (_seq->hsx.contigFilePos > hsxMaxFilePos)
		suicidef ("in load_hsx_sequence for %s,"
		          " file pos for contig %u (%010lX) exceeds max (%010lX)",
				  sequence_filename(_seq), _seq->contig,
				  _seq->hsx.contigFilePos, hsxMaxFilePos);
	err = fseek (_seq->f, (long int) _seq->hsx.contigFilePos, SEEK_SET);
	if (err != 0)
		suicidef ("in load_hsx_sequence for %s, index fseek(%010lX) returned %d",
		          sequence_filename(_seq), _seq->hsx.contigFilePos, err);

	seqName = read_hsx_index_entry (_seq);
	_seq->hsx.contigFilePos = (u64) ftell (_seq->f);

	if (!keeper) return;

	if (_seq->hsx.seqLength > (u64) maxSequenceLen)
		suicidef ("in load_hsx_sequence for %s, "
		          "sequence length " unsposFmt " for %s "
		          "exceeds maximum (" unsposFmt ")",
		          sequence_filename(_seq), _seq->hsx.seqLength, seqName,
		          maxSequenceLen);

	// copy the sequence name as our header

	numChars = strlen (seqName);
	if (_seq->headerSize < (unsigned) (numChars+1))
		{
		_seq->header      = realloc_or_die ("load_hsx_sequence (header)",
						 				  _seq->header, numChars+1);
		_seq->headerSize  = numChars+1;
		_seq->headerOwner = _seq->shortHeaderOwner = true;
		}

	strcpy (/*to*/ _seq->header, /*from*/ seqName);

	//////////
	// make sure we have enough room for the sequence's data
	// 
	// $$$ note that the allocated vector may be way more than needed, if start
	// $$$ .. and end limits reduce the number of characters we actually want;
	// $$$ .. this could be improved if it causes problems
	//////////

	dnaSize = (unspos) _seq->hsx.seqLength;

	if ((_seq->startLimit != 0) && (_seq->startLimit > dnaSize))
		suicidef ("beyond end in %s/%s (%ld > " unsposFmt ")",
		          sequence_filename(_seq), seqName, _seq->startLimit, dnaSize);

	if ((_seq->endLimit != 0) && (_seq->endLimit > dnaSize))
		{
		if (_seq->endIsSoft)
			{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
		else
			suicidef ("beyond end in %s/%s (%ld > " unsposFmt ")",
			          sequence_filename(_seq), seqName, _seq->endLimit, dnaSize);
		}

	startLimit = (u32) _seq->startLimit;
	endLimit   = (u32) _seq->endLimit;

	if (startLimit == 0) startLimit = 1;
	if (endLimit   == 0) endLimit   = dnaSize;

	_seq->trueLen += dnaSize;

#if (maxSequenceIndex <= 32)	// otherwise compiler complains that this test is
								// .. always false
	if ((dnaSize > maxSequenceLen) || (_seq->len > maxSequenceLen - dnaSize))
		suicidef ("in load_hsx_sequence for %s/%s, "
		          "sequence length %s+%s exceeds maximum (%s)",
		          sequence_filename(_seq), seqName,
		          commatize(_seq->len),commatize(dnaSize),
		          commatize(maxSequenceLen));
#endif

	sequence_long_enough (_seq, _seq->len + dnaSize, false);

	//////////
	// read the sequence's data
	//////////

	seqFName = _seq->hsx.fileInfo[_seq->hsx.seqFileIx].name;
	seqF     = _seq->hsx.fileInfo[_seq->hsx.seqFileIx].f;
	if (seqF == NULL)
		{
		// $$$ we should probably keep track of the number of open files and
		// $$$ close some (by LRU) if too many are open
		seqF = fopen_or_die (seqFName, "rb");
		_seq->hsx.fileInfo[_seq->hsx.seqFileIx].f = seqF;
		}

	if (_seq->hsx.seqFilePos > hsxMaxFilePos)
		suicidef ("in load_hsx_sequence for %s/%s,"
		          " file pos for sequence %s (%010lX) exceeds max (%010lX)",
				  sequence_filename(_seq), seqName, _seq->header,
				  _seq->hsx.seqFilePos, hsxMaxFilePos);
	err = fseek (seqF, _seq->hsx.seqFilePos, SEEK_SET);
	if (err != 0)
		suicidef ("in load_hsx_sequence for %s/s,"
		          " data fseek(%s,%08lX) returned %d",
		          sequence_filename(_seq), seqName,
		          seqFName, _seq->hsx.seqFilePos, err);

	// if the first character is a '>' (and the length is non-zero), we have to
	// skip this sequence header
	// $$$ it might be a good idea to validate that the header we are reading
	// $$$ .. matches the name of the sequence we think we're going to be reading

	prevCh = '\n';

	ch = getc_or_die (seqF, seqFName);
	if ((ch == '>') && (dnaSize != 0))
		{
		while (ch != '\n') // (skip line)
			ch = getc_or_die (seqF, seqFName);
		// get first character of next line
		ch = getc_or_die (seqF, seqFName);
		}

	while ((ch == ' ') || (ch == '\t')) // (skip whitespace)
		ch = getc_or_die (seqF, seqFName);

	// scan the file, keeping characters that are (a) nucleotides and (b) are
	// within our index limits

	index = 0;
	while (ch != EOF)
		{
		if ((prevCh == '\n') && (ch == '>')) // (start of next sequence)
			break;

		switch (char_to_fasta_type[(u8)ch])
			{
			case _nucleotide:
				break;
			case _ambiguous:
				if (!_seq->allowAmbiDNA) goto bad_char;
				break;
			case _newline:
				ch = '\n'; // (allow for unix, mac, or pc line ends)
				goto next_char;
			case _bad:
				goto bad_char;
			}

		// this is a nucleotide, do we want it?

		index++;

		if ((startLimit != 0) && (index < startLimit)) goto next_char;
		if ((endLimit   != 0) && (index > endLimit))   goto next_char;

		// ok, let's store it

		_seq->v[_seq->len++] = ch;

		// go try the next character

	next_char:
		prevCh = ch;
		do // (skip whitespace)
			{
			ch = getc_or_die (seqF, seqFName);
			} while ((ch == ' ') || (ch == '\t'));
		}

	_seq->v[_seq->len] = 0;				// (set the terminating zero)

	// $$$ we should make sure the sequence was as long as it said it was

	_seq->hsx.contigLoaded = true;

	if (startLimit == 0) _seq->start = 1;
	                else _seq->start = startLimit;

	return;

	// failure exits
	// $$$ report line number and sequence name here

bad_char:
	if (dna_isprint(ch))
		suicidef ("bad fasta character in %s, %s: %c",
		          sequence_filename(_seq), seqName, (int) ch);
	else
		suicidef ("bad fasta character in %s, %s (ascii %02X)",
		          sequence_filename(_seq), seqName, (u8) ch);
	}

//--- lookup_hsx_sequence ---

static u64 lookup_hsx_sequence
   (seq*	_seq,
	char*	name)
	{
	u32		bucket;
	u64		fileOffset;
	u64		bucketStart, bucketEnd;
	int		err;

	bucket = hassock_hash (name, strlen(name)) % _seq->hsx.numBuckets;
	fileOffset = _seq->hsx.hashTableOffset + (5 * (u64) bucket);
	if (fileOffset > hsxMaxFilePos)
		suicidef ("in lookup_hsx_sequence for %s,"
				  " file pos for %s hash bucket %d (%010lX) exceeds max (%010lX)",
				  sequence_filename(_seq), bucket, fileOffset, hsxMaxFilePos);
	err = fseek (_seq->f, (long int) fileOffset, SEEK_SET);
	if (err != 0)
		suicidef ("in lookup_hsx_sequence for %s, file table fseek(%010lX) returned %d",
		          sequence_filename(_seq), fileOffset, err);

	bucketStart = read_5 (_seq, _seq->hsx.bigEndian);
	if ((bucketStart & hsxMsBit5) != 0) // (bucket is empty)
		return hsxMsBit5; // (not found)
	if (bucketStart > hsxMaxFilePos)
		suicidef ("in lookup_hsx_sequence for %s,"
				  " file pos for %s bucket start (%010lX) exceeds max (%010lX)",
				  sequence_filename(_seq), bucketStart, hsxMaxFilePos);

	bucketEnd = read_5 (_seq, _seq->hsx.bigEndian) & ~hsxMsBit5;
	if (bucketEnd > hsxMaxFilePos)
		suicidef ("in lookup_hsx_sequence for %s,"
				  " file pos for %s bucket end (%010lX) exceeds max (%010lX)",
				  sequence_filename(_seq), bucketEnd, hsxMaxFilePos);

	return find_hsx_sequence (_seq, name, bucketStart, bucketEnd);
	}

//--- find_hsx_sequence ---

static u64 find_hsx_sequence
   (seq*	_seq,
	char*	name,
	u64		bucketStart,
	u64		bucketEnd)
	{
	u64		bucketOffset = bucketStart;
	char*	seqName;
	int		diff, err;

	err = fseek (_seq->f, (unsigned long) bucketOffset, SEEK_SET);
	if (err != 0)
		suicidef ("in find_hsx_sequence for %s,"
		          " file table fseek(%010lX) returned %d",
		          sequence_filename(_seq), bucketOffset, err);

	while (bucketOffset < bucketEnd)
		{
		seqName = read_hsx_index_entry (_seq);
		diff = strcmp (seqName, name);
		if (diff == 0) return bucketOffset; // (sequence name found)
		if (diff >  0) break;               // (sequence name not found)
		bucketOffset += 1 + 6 + 5 + strlen(seqName) + 1;
		}

	return hsxMsBit5; // (not found)
	}

//--- read_hsx_index_entry ---

static char* read_hsx_index_entry
   (seq* _seq)
	{
	_seq->hsx.seqLength  = read_5   (_seq, _seq->hsx.bigEndian);
	_seq->hsx.seqFileIx  = seq_getc (_seq);
	_seq->hsx.seqFilePos = read_6   (_seq, _seq->hsx.bigEndian);
	return read_hsx_string (_seq, _seq->f);
	}

//--- read_hsx_string ---

static char* read_hsx_string
   (seq*			_seq,
	FILE*			f)
	{
	static char		s[256];
	unsigned int	stringSize;
	size_t			bytesRead;

	// read the name

	stringSize = getc_or_die (_seq->f, _seq->fileName);
	if (stringSize == 0)
		{ s[0] = 0;  return s; }

	bytesRead = fread (s, 1, stringSize, f);
	if (bytesRead != stringSize)
		suicidef ("in read_hsx_string for %s, short read\n"
				  "wanted %d bytes, only got %d",
				  sequence_filename(_seq), stringSize, bytesRead);

	s[stringSize] = 0;
	return s;
	}

//----------
//
// load_qdna_sequence--
//	Load a quantum-dna sequence from the associated file.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//	int		keeper:	true  => actually load the sequence
//					false => just skip the sequence
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------
//
// Qdna file format:
//
//	Fields can be in big- or little-endian format;  they must match the
//	endianess of the magic number.
//
//	Version 2 (name ignored; files with named properites not supported):
//
//	offset 0x00: C4 B4 71 97   big endian magic number (97 71 B4 C4 => little endian)
//	offset 0x04: 00 00 02 00   version 2.0 (fourth byte is sub version)
//	offset 0x08: 00 00 00 14   header length (in bytes, including this field)
//	offset 0x0C: xx xx xx xx   S, offset (from file start) to data sequence
//	offset 0x10: xx xx xx xx   N, offset to name, 0 indicates no name
//	offset 0x14: xx xx xx xx   length of data sequence (counted in 'items')
//	offset 0x18: xx xx xx xx   (for version >= 2.0) P, offset to named
//							   .. properties, 0 indicates no properties
//	offset    N: ...           name (zero-terminated string)
//	offset    S: ...           data sequence
//	offset    P: ...           named properties (see below)
//
//	The named properties section is not allowed in this implementation.
//
//	Version 1 (name ignored):
//
//	offset 0x00: C4 B4 71 97   big endian magic number (97 71 B4 C4 => little endian)
//	offset 0x04: 00 00 01 00   version (fourth byte will be sub version)
//	offset 0x08: 00 00 00 10   header length (in bytes, including this field)
//	offset 0x0C: xx xx xx xx   S, offset (from file start) to data sequence
//	offset 0x10: xx xx xx xx   N, offset to name, 0 indicates no name
//	offset 0x14: xx xx xx xx   length of data sequence (counted in 'items')
//	offset    N:  ...          name (zero-terminated string)
//	offset    S:  ...          data sequence
//
//	Version 0:
//
//	offset 0x00: 9E 65 56 F6   magic number
//	offset 0x04:  ...          data sequence
//
//	Additionally, we will accept any binary file and interpret it as the data
//	sequence.  Note that if the data sequence happens to begin with one of the
//	magic numbers above, we will fail to read the file properly.  Further, if
//	the file contains newlines that are not part of the sequence, we will fail
//	to read the file properly.
//
//----------

// $$$ why don't we use the name from the file?????

static void load_qdna_sequence
   (seq*	_seq,
	int		keeper)
   	{
	u32		magic, version, headerLen, seqOffset, nameOffset, propOffset;
	u32		length, startLimit, startIndex, endLimit;
	unspos	newSeqLen;
	int		oldFormat, bigEndian, lengthKnown;
	int		ch;
	int		err, numChars;

	if (!keeper) return;
	if (_seq == NULL) suicide ("load_qdna_sequence(NULL)");

	//////////
	// process the header
	//////////

	// validate the magic number

	oldFormat = bigEndian = false;

	magic = read_4_big (_seq);
	if      (magic == qdnaMagicLittle)    { ; }
	else if (magic == qdnaMagicBig)       { bigEndian = true; }
	else if (magic == oldQdnaMagicLittle) { oldFormat = true; }
	else if (magic == oldQdnaMagicBig)    { oldFormat = bigEndian = true; }
	else
		{
		seq_ungetc ((magic >> 24) & 0xFF, _seq);
		seq_ungetc ((magic >> 16) & 0xFF, _seq);
		seq_ungetc ((magic >> 8)  & 0xFF, _seq);
		seq_ungetc ( magic        & 0xFF, _seq);
		oldFormat = true;
		}

	// skip the header (unless it's the old format)

	if (oldFormat)
		{
		lengthKnown = false;
		length      = 0;
		}
	else
		{
		version = read_4 (_seq, bigEndian);
		if (((version >> 8) != 1) && ((version >> 8) != 2))
			suicidef ("unsupported qdna version in %s (%08lX)",
			          sequence_filename(_seq), version);

		headerLen  = read_4 (_seq, bigEndian);
		seqOffset  = read_4 (_seq, bigEndian);
		nameOffset = read_4 (_seq, bigEndian);
		length     = read_4 (_seq, bigEndian);  lengthKnown = true;

		if ((version >> 8) == 1)
			skip_chars (_seq, seqOffset - 0x18);
		if ((version >> 8) == 2)
			{
			propOffset = read_4 (_seq, bigEndian);
			if (propOffset != 0)
				suicidef ("qdna named properties are not supported in %s",
				          sequence_filename(_seq));
			skip_chars (_seq, seqOffset - 0x1C);
			}

		_seq->trueLen += length;
		}

	//////////
	// skip ahead to the first desired base, and try to determine how many
	// bases we'll read
	//////////

	if (lengthKnown)
		{
		if ((_seq->startLimit != 0) && (_seq->startLimit > (unspos) length))
			suicidef ("beyond end in %s (%ld > %ld)",
			          sequence_filename(_seq), _seq->startLimit, length);

		if ((_seq->endLimit != 0) && (_seq->endLimit > (unspos) length))
			{
			if (_seq->endIsSoft)
				{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
			else
				suicidef ("beyond end in %s (%ld > %ld)",
				          sequence_filename(_seq), _seq->endLimit, length);
			}
		}
	else
		{
		if ((_seq->startLimit != 0) && (_seq->startLimit > (unspos) 0xFFFFFFFF))
			suicidef ("invalid start limit in %s (%ld > %ld)",
			          sequence_filename(_seq), _seq->startLimit, 0xFFFFFFFF);

		if ((_seq->endLimit != 0) && (_seq->endLimit > (unspos) 0xFFFFFFFF))
			{
			if (_seq->endIsSoft)
				{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
			else
				suicidef ("invalid end limit in %s (%ld > %ld)",
				          sequence_filename(_seq), _seq->endLimit, 0xFFFFFFFF);
			}
		}

	startLimit = (u32) _seq->startLimit;
	endLimit   = (u32) _seq->endLimit;

	if (startLimit == 0) startLimit = 1;
	startIndex = startLimit - 1;

	if (startIndex > 0)
		{
		if (!skip_chars (_seq, startIndex))
			suicidef ("bad start index for %s: %d",
			          sequence_filename(_seq), startIndex);
		}

	if (endLimit != 0)
		{ length = endLimit - startIndex;  lengthKnown = true; }

	//////////
	// allocate the vector (if we know the length)
	//////////

	newSeqLen = 0;
	if (lengthKnown)
		{
#if (maxSequenceIndex <= 32)	// otherwise compiler complains that this test is
								// .. always false
		if ((length > maxSequenceLen) || (_seq->len > maxSequenceLen - length))
			suicidef ("in load_qdna_sequence for %s, "
			          "sequence length %s+%s exceeds maximum (%s)",
			          sequence_filename(_seq),
			          commatize(_seq->len), commatize(length),
			          commatize(maxSequenceLen));
#endif

		newSeqLen = _seq->len + length;
		sequence_long_enough (_seq, newSeqLen, false);
		}

	//////////
	// read the sequence
	//////////

	while (true)
		{
		if ((newSeqLen != 0) && (_seq->len >= newSeqLen))
			break;

		// read the next character from the sequence

		ch = seq_getc (_seq);
		if (ch == EOF) break;

		if (ch == 0)
			suicidef ("in load_qdna_sequence(), file contains a zero");

		// allocate more room in the vector if we need it, and deposit the
		// character in the sequence

		if (!lengthKnown)
			sequence_long_enough (_seq, _seq->len+1, true);

		_seq->v[_seq->len++] = (u8) ch;
		}

	_seq->v[_seq->len] = 0;

	if (oldFormat)
		_seq->trueLen += _seq->len + startIndex;
	else
		{ ; } // (for new format, we already added the file length to _seq->trueLen

	if ((newSeqLen != 0) && (_seq->len < newSeqLen))
		suicidef ("beyond end in %s (%ld > end of file)",
		          sequence_filename(_seq), endLimit);

	_seq->pendingLen   = 0;					// (discard any pending chars)
	_seq->pendingStack = _seq->pendingChars + seqBufferSize;

	if (startLimit == 0) _seq->start = 1;
	                else _seq->start = startLimit;

	// skip the rest of the sequence

	if ((oldFormat) && (_seq->needTrueLen))
		{
		while ((ch = seq_getc (_seq)) != EOF)
			_seq->trueLen++;
		}
	else
		{
		err = fseek (_seq->f, 0, SEEK_END);	// skip the rest of the sequence
		if (err != 0)
			suicidef_with_perror ("in load_qdna_sequence(), fseek returned %d",
			                      err);
		}

	//////////
	// create a sequence header
	//////////

	if (!_seq->lockedHeader)
		{
		numChars = snprintf (_seq->header, 0, "%s:" unsposDashFmt,
		                     sequence_filename(_seq),
		                     _seq->start, _seq->start + _seq->len-1);

		if (_seq->headerSize < (unsigned) numChars+1)
			{
			_seq->header      = realloc_or_die ("load_qdna_sequence (header)",
			                                   _seq->header, numChars+1);
			_seq->headerSize  = numChars+1;
			_seq->headerOwner = _seq->shortHeaderOwner = true;
			}

		snprintf (_seq->header, numChars+1, "%s:" unsposDashFmt,
		          sequence_filename(_seq),
		          _seq->start, _seq->start + _seq->len-1);
		}
	}

//----------
//
// another_sequence--
//	Determine if the associated file has another sequence.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to dispose of.
//
// Returns:
//	true if there is another sequence to load, false if not.
//
//----------

int another_sequence
   (seq*			_seq)
	{
	debugNamesFile_3;

	if (_seq == NULL)                      suicide ("another_sequence(NULL)");
	if (_seq->fileType == seq_type_nofile) return false;

	// if we're subsampling the file's sequences, skip past sequences as needed

	if (_seq->subsampleN > 0)
		{
		if (_seq->subsampleSkip > 0)
			skip_sequences (_seq, _seq->subsampleSkip);
		_seq->subsampleSkip = 0;
		}

	return another_sequence_core (_seq);
	}

static int another_sequence_core
   (seq*			_seq)
	{
	seqpartition*	sp;
	int				ch;

	// if this is a partitioned sequence and we've finished loading it, there's
	// never another sequence

	sp = &_seq->partition;
	if ((sp->p != NULL) && (sp->state >= seqpart_ready))
		return false;

	// if we've previously positioned the file to a contig, but have yet to
	// load that contig, then we have another sequence

	if (_seq->contigPending) return true;

	// if we have a contigs-of-interest file (and we haven't pre-loaded a
	// sequence), get the next contig name

	if ((_seq->namesFile != NULL) && (!_seq->preLoaded))
		{
		if (!read_contig_name (_seq))
			return false;
		}

	// for 2bit or hsx files it's a matter of whether we're at the end of the
	// index list

	if (_seq->fileType == seq_type_2bit)
		{
		if (!_seq->twoBit.contigLoaded)              return (_seq->twoBit.numContigs > 0);
		if (_seq->contigOfInterest != NULL)          return false;
		if (_seq->namesFile != NULL)                 return find_next_2bit_coi (_seq);
		if (_seq->contig >= _seq->twoBit.numContigs) return false;
		return true;
		}
	else if (_seq->fileType == seq_type_hsx)
		{
		if (!_seq->hsx.contigLoaded)                 return (_seq->hsx.numContigs > 0);
		if (_seq->contigOfInterest != NULL)          return false;
		if (_seq->namesFile != NULL)                 return find_next_hsx_coi (_seq);
		if (_seq->contig >= _seq->hsx.numContigs)    return false;
		return true;
		}

	// if we pre-loaded a sequence then that sequence counts as "another"
	// sequence since the caller doesn't know we did so
	// $$$ why isn't this check *above* the 2bit/hsx tests?

	if (_seq->preLoaded) return true;

	// otherwise it's a matter of having data left in the file

	if (feof   (_seq->f)) return false;		// we've previously hit end of file
	if (ferror (_seq->f)) return false;		// we've previously had a problem

	if ((_seq->fileType == seq_type_fasta)	// we've got the next contig-of-
	 && (_seq->namesFile != NULL))			// .. interest
		return find_next_fasta_coi (_seq);
	else if ((_seq->fileType == seq_type_csfasta)
	 && (_seq->namesFile != NULL))
		return find_next_csfasta_coi (_seq);

	if (_seq->pendingLen > 0) return true;	// we have characters to process

	ch = getc_or_die (_seq->f,				// take a peek and see what's left
	                  _seq->fileName);
	if (ch == EOF) return false;			// we're at end of file now

	seq_ungetc (ch, _seq);					// save what we peeked at
	return true;							// we have characters to process
	}


// find_next_fasta_coi, find_next_csfasta_coi--
//	advance to the next contig-of-interest in fasta or csfasta file
//  (always returns true)

static int find_next_fasta_coi (seq* _seq)
	{ return find_next_general_fasta_coi (_seq, false); }

static int find_next_csfasta_coi (seq* _seq)
	{ return find_next_general_fasta_coi (_seq, true); }

static int find_next_general_fasta_coi
   (seq*	_seq,
	int		allowComments)
	{
	char	buffer[maxSequenceHeader+1];
	char*	header;
	int		headerLen;
	int		mustBeHeader;
	int		leadingWhite;
	char	ch, *s;
	int		ix;

	debugNamesFile_4;

	mustBeHeader = true;

	while (true)
		{
		// find the next header

		ch = seq_getc (_seq);
		if (ch == EOF) goto failure;

		if ((allowComments) && (ch == '#'))
			{ // comment, skip to end-of-line and go back and try again
			while (ch != '\n')
				{
				ch = seq_getc (_seq);
				if (ch == EOF) goto failure;
				}
			continue;
			}

		if (ch != '>')
			{
			if (mustBeHeader)
				suicidef ("internal error in find_next_fasta_coi\n"
				          "processing %s, looking for \"%s\"\n",
				          sequence_filename(_seq), _seq->nextContigName);
			continue;
			}

		if (!mustBeHeader) _seq->contig++;
		mustBeHeader = false;

		// skip leading white space

		debugNamesFile_5;

		leadingWhite = 0;

		ch = seq_getc (_seq);
		if (ch == EOF) goto failure;
		while ((ch != '\n') && (isspace (ch)))
			{
			leadingWhite++;
			ch = seq_getc (_seq);
			if (ch == EOF) goto failure;
			}

		if (ch == '\n')
			continue;	// (unnamed sequence)

		// read the header

		s = buffer;
		while (ch != '\n')
			{
			if (s - buffer >= maxSequenceHeader)	// (overflow;
				break;								//  .. truncate the header)
			*(s++) = ch;
			ch = seq_getc (_seq);
			if (ch == EOF) goto failure;
			}
		*s = 0;

		// if we have a name trigger, locate the sequence's name

		if (_seq->nameTrigger != NULL)
			{
			header = strstr (buffer, _seq->nameTrigger);
			if (header == NULL) continue;	// (effectively unnamed sequence)
			header += strlen (_seq->nameTrigger);

			s = header;
			while ((*s != 0) && ((isalnum(*s)) || (*s == '_')))
				s++;
			headerLen = s-header;
			}
		else if (!_seq->useFullNames)
			{
			shorten_header (/* from */ buffer, _seq->nameParseType, false,
			                /* to   */ NULL, NULL);
			header    = buffer;
			headerLen = strlen(buffer);
			}
		else
			{
			header    = buffer;
			headerLen = strlen(buffer);
			}

		// compare header to the contig-of-interest

		if (strncmp (header, _seq->nextContigName, headerLen) != 0)
			continue;
		if ((int) strlen (_seq->nextContigName) != headerLen)
			continue;

		break; // found a match!
		}

	debugNamesFile_6;

	// unget the header

	seq_ungetc (ch, _seq);	// (ch terminated the header)

	for (ix=strlen(buffer) ; ix>0 ; )
		seq_ungetc (buffer[--ix], _seq);

	while (leadingWhite-- > 0) seq_ungetc (' ', _seq);
	seq_ungetc ('>', _seq);

	_seq->contigPending = true;
	return true;

	// failure, the contig name was not found

failure:
	suicidef ("%s does not contain (or contains out of order)\n"
	          "         the contig \"%s\"",
			  sequence_filename(_seq), _seq->nextContigName);
	return false; // (will never reach here)
	}


// find_next_2bit_coi--
//	advance to the next contig-of-interest in 2bit header
//  (always returns true)

static int find_next_2bit_coi
   (seq*		_seq)
	{
	char		seqName[maxSequenceName+1];
	long int	savedContigFilePos, seqDataPos;
	int			err;

	debugNamesFile_7;

	// position to the sequence's next index table entry

	err = fseek (_seq->f, _seq->twoBit.contigFilePos, SEEK_SET);
	if (err != 0)
		suicidef ("in find_next_2bit_coi(%s), index fseek(%08lX) returned %d",
		          sequence_filename(_seq), _seq->twoBit.contigFilePos, err);

	// read index table entries until we find the one we're looking for

	while (true)
		{
		if (_seq->contig >= _seq->twoBit.numContigs)
			suicidef ("%s does not contain (or contains out of order)\n"
			          "         the contig \"%s\"",
			          sequence_filename(_seq), _seq->nextContigName);

		// read the sequence's next index table entry

		savedContigFilePos = ftell (_seq->f);
		seqDataPos = read_2bit_index_entry (_seq, seqName, _seq->contig);
		if (strcmp (seqName, _seq->nextContigName) == 0) break;

		_seq->contig++;
		}

	_seq->twoBit.contigFilePos = savedContigFilePos;
	_seq->contigPending        = true;
	return true;
	}


// find_next_hsx_coi--
//	advance to the next contig-of-interest in hsx header
//  (always returns true)

static int find_next_hsx_coi
   (seq*	_seq)
	{
	u64		fileOffset;

	debugNamesFile_8;

	fileOffset = lookup_hsx_sequence (_seq, _seq->nextContigName);

	if ((fileOffset & hsxMsBit5) != 0)
		suicidef ("hsx file %s doesn't contain %s",
				  sequence_filename(_seq), _seq->nextContigName);
	if (fileOffset > hsxMaxFilePos)
		suicidef ("in find_next_hsx_coi for %s,"
				  " file pos for %s (%010lX) exceeds max (%010lX)",
				  sequence_filename(_seq), _seq->nextContigName, fileOffset);

	_seq->hsx.contigFilePos = fileOffset;
	_seq->contigPending     = true;
	return true;
	}

//----------
//
// read_contig_name--
//	Read the next name from a contigs-of-interest file.
//
// The file will contain one name per line.  Any leading whitespace is ignored,
// any comment lines are ignored (# is the comment character), and the name is
// only up to the first whitespace character.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence.
//
// Returns:
//	true if we were successful;  false if there are no more names in the file.
//
//----------

static int read_contig_name
   (seq*	_seq)
	{
	char*	line     = _seq->nextContigName;
	int		lineSize = sizeof(_seq->nextContigName);
	char	discard[maxSequenceName+1];
	int		len;
	int		missingEol;
	char*	waffle, *s;

	while (fgets (line, lineSize, _seq->namesFile) != NULL)
		{
		// check for lines getting split by fgets (the final line in the file
		// might not have a newline, but no internal lines can be that way);
		// if the line was split we simply read ahead until we find the end of
		// the line

		len = strlen(line);
		if (len == 0) continue;
		missingEol = (line[len-1] != '\n');

		if (missingEol)
			{
			while (fgets (discard, sizeof(discard), _seq->namesFile) != NULL)
				{
				len = strlen(discard);
				if (len == 0) break;
				if (discard[len-1] == '\n') break;
				}
			}

		// trim blanks, end of line, and comments, and ignore blank lines

		len = strlen(line);
		if (line[len-1] == '\n') line[--len] = 0;

		waffle = strchr (line, '#');
		if (waffle != NULL) *waffle = 0;

		trim_string (line);
		if (line[0] == 0) continue;

		// ok, the line has something in it

		s = skip_darkspace(line);
		*s = 0;

		return true;
		}

	return false;
	}

//----------
//
// create_short_header--
//	Convert a sequence's header into a shorter version.  The shorter version is
//	intended to be useful as a sequence's name in maf or axt files.
//
// Examples:
//
//	>/depot/data1/cache/human/hg18/_seq/chr16.nib:120000-190000   chr16
//  owl_monkey 122000-180000 of owl_monkey.ENm008.fa             owl_monkey
//  > armadillo|ENm001|JAN-2006|9361|NISC|...|1|1|.              armadillo
//	>reverse complement of /depot/../chr14.nib:...               chr14
//	>positions 180000-250000 of armadillo|ENm008|...             armadillo
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence.
//
// Returns:
//	(nothing;  _seq->shortHeader and _seq->shortHeaderSize are modified)
//
//----------
//
// Note:	It is possible for the resulting short header name to be an empty
//			string.
//
//----------

static void create_short_header
   (seq*	_seq)
	{
	int		skipPath;

	if (_seq->header == NULL)
		{
		if ((_seq->shortHeader != NULL) && (_seq->shortHeaderSize != 0))
			_seq->shortHeader[0] = 0;
		return;
		}

	if ((!_seq->headerOwner) || (!_seq->shortHeaderOwner))
		{
		char* name = (_seq->fileName != NULL)? _seq->fileName
		                                     : _seq->header;
		suicidef ("internal error, attempt to shorten external sequence header (%s)",
		          name);
		}

	if (strstr(_seq->header,"{number}") != NULL)
		expand_nickname( /* from */ _seq->header, _seq->contig,
	                     /* to   */ &_seq->shortHeader, &_seq->shortHeaderSize);
	else
		{
		skipPath = (_seq->fileType == seq_type_nib);
		shorten_header (/* from */ _seq->header, _seq->nameParseType, skipPath,
		                /* to   */ &_seq->shortHeader, &_seq->shortHeaderSize);
		}
	}


static void shorten_header
   (char*	src,
	int		nameParseType,
	int		skipPath,
	char**	_dst, // (NULL => write it in place)
	u32*	_dstSize)
	{
	char*	dst;
	u32		dstSize;
	char*	h, *hh, *s;
	u32		len, sLen;

	// skip fasta '>', leading whitespace, and/or a path

	h = src;
	if (h[0] == '>') h++;
	h = skip_whitespace (h);

	// skip "reverse complement of" and/or "positions A-B of"

	s = "reverse complement of ";
	if (strcmp_prefix (h, s) == 0)
		h = skip_whitespace (h + strlen (s));

	s = "positions ";
	if (strcmp_prefix (h, s) == 0)
		{
		hh = skip_whitespace (h + strlen (s));
		hh = skip_darkspace (hh);
		hh = skip_whitespace (hh);
		s = "of ";
		if (strcmp_prefix (hh, s) == 0)				// (we only change h if
			h = skip_whitespace (hh + strlen (s));	//  .. "of" is present)
		}

	// skip a path

	if (skipPath)
		{
		while (true)
			{
			hh = strchr (h, pathSlash);
			if (hh == NULL) break;
			h = hh + 1;
			}
		}

	h = skip_whitespace (h);

	// figure out the length to copy;  we'll truncate at the first space or
	// "funny" character

	if (nameParseType == name_parse_type_alnum)
		{
		len = strspn (h, "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		                 "abcdefghijklmnopqrstuvwxyz"
		                 "0123456789" "_");
		goto skip_suffix_trim;
		}
	else if (nameParseType == name_parse_type_darkspace)
		len = strcspn (h, " \t");
	else // if (nameParseType == name_parse_type_core)
		len = strcspn (h, " \t|:");

	// if the suffix is ".nib", ".fasta", etc., remove it

	s = ".nib";
	sLen = strlen(s);
	if ((len > sLen) && (strncmp (h+len-sLen,s,sLen) == 0))
		len -= sLen;

	s = ".2bit";
	sLen = strlen(s);
	if ((len > sLen) && (strncmp (h+len-sLen,s,sLen) == 0))
		len -= sLen;

	s = ".hsx";
	sLen = strlen(s);
	if ((len > sLen) && (strncmp (h+len-sLen,s,sLen) == 0))
		len -= sLen;

	s = ".fasta";
	sLen = strlen(s);
	if ((len > sLen) && (strncmp (h+len-sLen,s,sLen) == 0))
		len -= sLen;

	s = ".fa";
	sLen = strlen(s);
	if ((len > sLen) && (strncmp (h+len-sLen,s,sLen) == 0))
		len -= sLen;

skip_suffix_trim:

	// create the header

	if (_dst == NULL)
		{
		strncpy (src, h, len);
		src[len] = 0;
		}
	else
		{
		dst     = *_dst;
		dstSize = *_dstSize;

		if (len+1 > dstSize)
			{
			dst = realloc_or_die ("shorten_header", dst, len+1);
			dstSize = len+1;
			}

		strncpy (dst, h, len);
		dst[len] = 0;
		*_dst     = dst;
		*_dstSize = dstSize;
		}
	}


static void expand_nickname
   (char*	src,
	u32		contigNumber,
	char**	_dst,
	u32*	_dstSize)
	{
	char*	dst;
	u32		dstSize;
	char*	s, *d, *expand;
	u32		len;

	// determine the size of the resulting header

	len = strlen (src)
	    - strlen ("{number}");
	    + snprintf (src, 0, unsposFmt, contigNumber);

	// allocate the header (if necessary)

	dst     = *_dst;
	dstSize = *_dstSize;

	if (len+1 > dstSize)
		{
		dst = realloc_or_die ("expand_nickname", dst, len+1);
		dstSize = len+1;
		}

	// create the header

	s = src;
	d = dst;

	expand = strstr (src, "{number}");
	if (expand > src)
		{
		strncpy (d, s, expand-src);
		d += expand-src;
		s =  expand + strlen("{number}");
		}

	sprintf (d, unsposFmt, contigNumber);
	d += strlen(d);

	strcpy (d, s);

	*_dst     = dst;
	*_dstSize = dstSize;
	}

//----------
//
// enough_partitions--
//	Make sure a sequence has enough room for a specified number of partitions.
//
//----------
//
// Arguments:
//	seq*	_seq:			The sequence to check.
//	u32		numPartitions:	The number of partitions to allocate for (not
//							.. including the extra partition used as a
//							.. sentinel).
//	u32		poolSize:		The number of bytes to allocate for a pool of
//							.. headers.  If this is zero, we will estimate the
//							.. pool size from the number of partitions.
//	int		anticipate:		true  => allocate extra, anticipating the need for
//							         .. more
//							false => don't
//
// Returns:
//	nothing;  _seq->partition.p and _seq->partition.pool may be modified;
//	failures result in fatality.
//
//----------

#define averageHeaderSize 20

static void enough_partitions
   (seq*			_seq,
	u32				numPartitions,
	u32				poolSize,
	int				anticipate)
	{
	seqpartition*	sp = &_seq->partition;
	u32				bytesNeeded;

	// if we have enough already, just return

	if ((sp->p != NULL)
	 && (sp->size > numPartitions)
	 && (sp->pool != NULL)
	 && ((poolSize > 0) && (sp->poolSize >= poolSize)))
		return;

	if (!sp->poolOwner)
		{
		char* name = (_seq->fileName != NULL)? _seq->fileName
		                                     : _seq->header;
		suicidef ("internal error, attempt to resize external partition names pool (%s)",
		          name);
		}

	if (sp->p    == NULL) sp->len     = 0;
	if (sp->pool == NULL) sp->poolLen = 0;

	if (poolSize == 0) poolSize = numPartitions * (averageHeaderSize + 1);

	// allocate partition array;  note that we bump up the number of records
	// allocated to as many as can fit in a multiple of 16K

	if (sp->size <= numPartitions)
		{
		numPartitions++;					// (extra one for a sentinel)
		if (anticipate)						// anticipatory, grow by about 13%
			numPartitions += 30 + numPartitions / 8;
		bytesNeeded   = numPartitions * sizeof(partition);
		bytesNeeded   = round_up_16K (bytesNeeded);
		numPartitions = bytesNeeded / sizeof(partition);
		bytesNeeded   = numPartitions * sizeof(partition);
		sp->p         = realloc_or_die ("enough_partitions (p)",
		                                sp->p, bytesNeeded);
		sp->size      = numPartitions;
		}

	// allocate pool for partition headers

	if (sp->poolSize < poolSize)
		{
		if (anticipate)						// anticipatory, grow by about 13%
			poolSize += 30*(averageHeaderSize+1) + poolSize / 8;
		bytesNeeded   = round_up_16K (poolSize);
		sp->pool      = realloc_or_die ("enough_partitions (pool)",
		                                sp->pool, bytesNeeded);
		sp->poolSize  = bytesNeeded;
		sp->poolOwner = true;
		}

	}

//----------
//
// lookup_partition--
//	Map a position in a partitioned sequence to its partition record.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence.
//	unspos	pos:	The position to look up.  This is an index into _seq->v,
//					.. (origin zero).
//
// Returns:
//	A pointer to the partition record (see note (1);  failures result in
//	fatality.
//
//----------
//
//	Notes:
//
//	(1)	The pointer p returned points to the partition record which contains the
//		contig, header, and sepPos for the partition.  sepPos is the bounding
//		NUL at the left end (lower index) of the sequence.  The pointer p+1
//		points to the next partition record, whose sepPos is the bounding NUL
//		at the right end (higher index) of the sequence. 
//
//----------

partition* lookup_partition
   (seq*			_seq,
	unspos			pos)
	{
	seqpartition*	sp = &_seq->partition;
	partition*		p;
	u32				hi, lo, ix;

	if (sp->p   == NULL) goto failure;
	if (sp->len == 0)    goto failure;

	p = sp->p;

	lo = 0;
	hi = sp->len;

	if (pos <= p[lo].sepPos) goto failure;
	if (pos >= p[hi].sepPos) goto failure;

	// perform binary search for the position
	//
	// loop invariants:					loop termination:
	//		0 <= lo < hi <= len				ix = lo = hi-1
	//		pos > p[lo].sepPos				pos > p[ix].sepPos
	//		pos < p[hi].sepPos				pos < p[ix+1].sepPos

	while (true)
		{
		ix = (lo + hi) / 2;			// when hi==lo+1, ix==lo
		if (hi == lo+1) break;
		if      (pos < p[ix].sepPos) hi = ix;
		else if (pos > p[ix].sepPos) lo = ix;
		else                         goto failure;	// pos == p[ix].sepPos,
													// .. which is illegal
		}

	// success

	return &sp->p[ix];

	// failure

failure:
	if (_seq->fileName == NULL)
		suicidef ("lookup_partition could not locate position " unsposFmt,
		          pos);
	else
		suicidef ("lookup_partition could not locate position " unsposFmt
		          " in %s",
		          pos, _seq->fileName);

	return NULL; // (never gets here)
	}

//----------
//
// lookup_named_partition--
//	Map a name to the correpsonding partition record in a partitioned sequence.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence.
//	name	name:	The name to look up.
//
// Returns:
//	A pointer to the partition record;  NULL if not found.
//
//----------

partition* lookup_named_partition
   (seq*			_seq,
	char*			name)
	{
	seqpartition*	sp = &_seq->partition;
	partition*		p;
	u32				ix;

	if (sp->p   == NULL) return NULL;
	if (sp->len == 0)    return NULL;

	// perform linear search for the name

	p = sp->p;

	for (ix=0 ; ix<sp->len ; ix++)
		{ if (strcmp (&sp->pool[p[ix].header], name) == 0) return &sp->p[ix]; }

	return NULL;
	}

//----------
//
// print_sequence--
//	Print a sequence to a file.
//
//----------
//
// Arguments:
//	FILE*	f:			The file to write to.
//	seq*	_seq:		The sequence to print.
//	char*	header:		The header to give the sequence.  If this is NULL,
//						.. _seq->header is used.  If this is empty, we assume
//						.. the caller has already taken care of printing the
//						.. header.
//	int		perLine:	The number of nucleotides to print per line.  Zero
//						.. indicates that they should all be printed on one
//						.. line.
//
// Returns:
//	(nothing)
//
//----------

void print_sequence
   (FILE*	f,
	seq*	_seq,
	char*	header,
	int		perLine)
	{
	unspos	ix;

	if (_seq == NULL)
		{ fprintf (f, "(null sequence)\n");  return; }

	if ((header == NULL) || (header[0] != 0))
		{
		if (header == NULL)
			{
			header = _seq->header;
			if (header == NULL)
				header = "";
			}

		if (header[0] == '>')
			header = skip_whitespace(header+1);

		if (header[0] == 0) fprintf (f, ">\n");
		               else fprintf (f, "> %s\n", header);
		}

	if (_seq->fileType == seq_type_qdna)
		{
		for (ix=0 ; ix<_seq->len ; ix++)
			{
			if ((ix != 0) && ((ix % perLine) == 0)) fprintf (f, "\n");
			fprintf (f, " %02X", _seq->v[ix]);
			}
		fprintf (f, "\n");
		}
	else
		{
		for (ix=0 ; ix<_seq->len ; ix++)
			{
			if ((ix != 0) && ((ix % perLine) == 0)) fprintf (f, "\n");
			fprintf (f, "%c", _seq->v[ix]);
			}
		fprintf (f, "\n");
		}
	}

//----------
//
// print_partition_table--
//	Print a sequence's partition table.  (for debugging only)
//
//----------
//
// Arguments:
//	FILE*	f:		The file to write to.
//	seq*	_seq:	The sequence to print.
//
// Returns:
//	(nothing)
//
//----------

#ifdef debugPartitions

static void print_partition_table
   (FILE*			f,
	seq*			_seq)
	{
	seqpartition*	sp = &_seq->partition;
	partition*		p;
	int				ix;

	sp = &_seq->partition;
	if (sp->p == NULL)
		{ fprintf (f, "sequence has no partition\n");  return; }
	if (sp->state == seqpart_empty)
		{ fprintf (f, "partition table is in empty state\n");  return; }
	if (sp->state == seqpart_loading)
		{ fprintf (f, "partition table is in loading state\n");  return; }

	p = sp->p;
	for (ix=0 ; ix<=sp->len ; ix++)
		{
		fprintf (f, "[%2d] %8u", ix, p[ix].sepPos);
		if (ix < sp->len)
			fprintf (f, " %3d %s", p[ix].contig, p[ix].header);
		fprintf (f, "\n");
		}

	}

#endif // debugPartitions

//----------
//
// mask_sequence, mask_sequence_keep--
//	Mask a sequence, in place, as prescribed by some file.
//
// mask_sequence() replaces any base in the prescribed intervals.
// mask_sequence_keep() replaces any base NOT in the prescribed intervals.
//
// A typical masking file looks like this:
//
//	1527933 3184039
//	4165389 6877343
//	7374477 7902860
//
// Each line describes a region to be masked.  Indexes are one-based, and
// inclusive on both ends.  Numbers are free format.  Comment lines (beginning
// with #) are ignored, as are blank lines.  Additional information after the
// first two columns is also ignored.
//
// Note that if the sequence has been reverse complemented, the masking
// intervals are relative to teh reverse strand.
//
//----------
//
// Arguments:
//	seq*	_seq:			The sequence to mask.
//	char*	maskFileName:	The name of the file containing masking
//							.. information.
//	int		maskChar:		The character to ask as a mask.  Normally this is a
//							.. character (e.g. 'X' or 'N').  However, the value
//							.. -1 means that we should mask by changing to
//							.. lowercase.
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

void mask_sequence
   (seq*	_seq,
	char*	maskFileName,
	int		_maskChar)
	{
	char	line[511+1], discard[511+1];
	char	maskChar = 0;
	int		toLower  = false;
	FILE*	maskF;
	int		len;
	int		lineNum, missingEol;
	char*	waffle;
	unspos	b, e;
	char	extra;
	int		numItems;

	if (_maskChar == -1) toLower  = true;
	else                 maskChar = (u8) _maskChar;

	// read the masking intervals and deposit the masking character at every
	// contained base

	maskF = fopen_or_die (maskFileName, "rt");

	lineNum = 0;
	while (fgets (line, sizeof(line), maskF) != NULL)
		{
		lineNum++;

		// check for lines getting split by fgets (the final line in the file
		// might not have a newline, but no internal lines can be that way);
		// if the line was split we simply read ahead until we find the end of
		// the line

		len = strlen(line);
		if (len == 0) continue;
		missingEol = (line[len-1] != '\n');

		if (missingEol)
			{
			while (fgets (discard, sizeof(discard), maskF) != NULL)
				{
				len = strlen(discard);
				if (len == 0) break;
				if (discard[len-1] == '\n') break;
				}
			}

		// trim blanks, end of line, and comments, and ignore blank lines

		len = strlen(line);
		if (line[len-1] == '\n') line[--len] = 0;

		waffle = strchr (line, '#');
		if (waffle != NULL) *waffle = 0;

		trim_string (line);
		if (line[0] == 0) continue;

		// ok, the line has something in it;  parse it as an interval

		numItems = sscanf (line, unsposFmtScanf " " unsposFmtScanf "%c",
		                   &b, &e, &extra);

		if ((numItems == 3) && ((extra == ' ') || (extra == '\t')))
			numItems = 2;
		if (numItems != 2)
			suicidef ("bad interval (in %s, line %d): \"%s\"\n",
			          maskFileName, lineNum, line);

		// trim the interval to our subsequence

		if (e < _seq->start) continue;
		if (b < _seq->start) b = _seq->start;
		b -= _seq->start;    // (nota bene, b is zero-based afterwards)
		e -= _seq->start-1;  // (nota bene, e is open interval afterwards)
		if (b >= _seq->len) continue;
		if (e >= _seq->len) e = _seq->len;

		// mask the interval

		if (e > b)
			{
			if (toLower)
				{
				for ( ; b<e ; b++)
					{
					if ((_seq->v[b] >= 'A') && (_seq->v[b] <= 'Z'))
						_seq->v[b] = _seq->v[b] + 'a' - 'A';
					}
				}
			else
				memset (_seq->v+b, maskChar, (size_t) (e-b));
			}
		}

	fclose(maskF);
	}

void mask_sequence_keep
   (seq*	_seq,
	char*	maskFileName,
	int		_maskChar)
	{
	char	line[511+1], discard[511+1];
	char	maskChar = 0;
	int		toLower  = false;
	FILE*	maskF;
	int		len;
	int		lineNum, missingEol;
	char*	waffle;
	unspos	b, e;
	char	extra;
	int		numItems;

	if ((_seq->fileType != seq_type_fasta)
	 && (_seq->fileType != seq_type_nib)
	 && (_seq->fileType != seq_type_2bit)
	 && (_seq->fileType != seq_type_hsx))
		suicidef ("masking of interval complements only valid for DNA sequences\n"
		          " (%s is %s file)",
		          sequence_filename(_seq), seqTypeNames[_seq->fileType]);

	if (_maskChar == -1) toLower  = true;
	else                 maskChar = (u8) _maskChar;

	// read the masking intervals and mark the most-significant bit of every
	// contained base

	maskF = fopen_or_die (maskFileName, "rt");

	lineNum = 0;
	while (fgets (line, sizeof(line), maskF) != NULL)
		{
		lineNum++;

		// check for lines getting split by fgets (the final line in the file
		// might not have a newline, but no internal lines can be that way);
		// if the line was split we simply read ahead until we find the end of
		// the line

		len = strlen(line);
		if (len == 0) continue;
		missingEol = (line[len-1] != '\n');

		if (missingEol)
			{
			while (fgets (discard, sizeof(discard), maskF) != NULL)
				{
				len = strlen(discard);
				if (len == 0) break;
				if (discard[len-1] == '\n') break;
				}
			}

		// trim blanks, end of line, and comments, and ignore blank lines

		len = strlen(line);
		if (line[len-1] == '\n') line[--len] = 0;

		waffle = strchr (line, '#');
		if (waffle != NULL) *waffle = 0;

		trim_string (line);
		if (line[0] == 0) continue;

		// ok, the line has something in it;  parse it as an interval

		numItems = sscanf (line, unsposFmtScanf " " unsposFmtScanf "%c",
		                   &b, &e, &extra);

		if ((numItems == 3) && ((extra == ' ') || (extra == '\t')))
			numItems = 2;
		if (numItems != 2)
			suicidef ("bad interval (in %s, line %d): \"%s\"\n",
			          maskFileName, lineNum, line);

		// trim the interval to our subsequence

		if (e < _seq->start) continue;
		if (b < _seq->start) b = _seq->start;
		b -= _seq->start;    // (nota bene, b is zero-based afterwards)
		e -= _seq->start-1;  // (nota bene, e is open interval afterwards)
		if (b >= _seq->len) continue;
		if (e >= _seq->len) e = _seq->len;

		// mark the interval

		while (b<e) _seq->v[b++] |= 0x80;
		}

	fclose(maskF);

	// scan the sequence, replacing unmarked bases with the masking character,
	// and removing the marks

	for (b=0 ; b<_seq->len ; b++)
		{
		if (_seq->v[b] == 0) continue;
		if ((_seq->v[b] & 0x80) != 0)	// marked => erase mark
			_seq->v[b] &= ~0x80;
		else if (!toLower)				// unmarked => mask it
			_seq->v[b] = maskChar;
		else							// unmarked => change to lower case
			{
			if ((_seq->v[b] >= 'A') && (_seq->v[b] <= 'Z'))
				_seq->v[b] = _seq->v[b] + 'a' - 'A';
			}
		}

	}

//----------
//
// colorize_sequence--
//	Create a sequence's color-space equivalent (see description below).
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to modify.
//
// Returns:
//	(nothing)
//
//----------
//
// We colorize by mapping each pair of nucleotides to a single "color" value
// according to the table below.  This extends the normal color space
// definition to account for lower case, N, and other characters, in a way that
// makes sense in the context of lastz's alignment operations.  The DNA sequence
// may contain upper and lower case nucleotides, N's, and other stuff.  Further,
// a paritioned sequence will contain NUL character separators (shown as '*' in
// the table below).  We also treat the DNA sequence as if it had an 'X'
// prepended to it.
//
//		                  second in pair
//		      | A | C | G | T | a | c | g | t | N/n | other | * |
//		------+---+---+---+---+---+---+---+---+-----+-------+---+
//		   A  | 0 | 1 | 2 | 3 | 0 | 1 | 2 | 3 |  N  |   X   | * |
//		   C  | 1 | 0 | 3 | 2 | 1 | 0 | 3 | 2 |  N  |   X   | * |
//		   G  | 2 | 3 | 0 | 1 | 2 | 3 | 0 | 1 |  N  |   X   | * |
//	first  T  | 3 | 1 | 2 | 0 | 3 | 1 | 2 | 0 |  N  |   X   | * |
//	 in	   a  | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |  N  |   X   | * |
//	  pair c  | 1 | 0 | 3 | 2 | 5 | 4 | 7 | 6 |  N  |   X   | * |
//		   g  | 2 | 3 | 0 | 1 | 6 | 7 | 4 | 5 |  N  |   X   | * |
//		   t  | 3 | 1 | 2 | 0 | 7 | 6 | 5 | 4 |  N  |   X   | * |
//		  N/n | N | N | N | N | N | N | N | N |  N  |   X   | * |
//		------+---+---+---+---+---+---+---+---+-----+-------+---+
//		other | X | X | X | X | X | X | X | X |  X  |   X   | * |
//		------+---+---+---+---+---+---+---+---+-----+-------+---+
//		   *  | X | X | X | X | X | X | X | X |  X  |   X   | * |
//		------+---+---+---+---+---+---+---+---+-----+-------+---+
//
// Example:
//
//	 dna:   G T C G A A C C C G * C A A C C G T A T T * T A A T A A G T T T
//	 color: X 1 2 3 2 0 1 0 0 3 * X 1 0 1 0 3 1 3 3 0 * X 3 0 3 3 0 2 1 0 0 
//
//----------

static void do_colorize (u8* seq, u8* colorSeq, unspos seqLen);

void colorize_sequence
   (seq*	_seq)
	{
	char*	name;

	if (_seq == NULL) suicide ("colorize_sequence(NULL)");
	if (_seq->len < 1) return;

	if (!_seq->vcOwner)
		{
		name = (_seq->fileName != NULL)? _seq->fileName : _seq->header;
		suicidef ("internal error, attempt to colorize external sequence (%s)",
		          name);
		}

	if (_seq->vc != NULL)
		{
		name = (_seq->fileName != NULL)? _seq->fileName : _seq->header;
		suicidef ("internal error, attempt to re-colorize sequence (%s)",
		          name);
		}

	_seq->vc = malloc_or_die ("colorize_sequence (vc)", _seq->len+1);
	do_colorize (_seq->v, _seq->vc, _seq->len);
	}

#define A_ 0
#define C_ 1
#define G_ 2
#define T_ 3
#define a_ 4
#define c_ 5
#define g_ 6
#define t_ 7
#define N_ 8
#define X_ 9
#define __ X_
#define Z_ 10

const s8 nuc_to_color_bits[256] =
	{
	Z_,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 0x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 1x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 2x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 3x (numbers)
	__,A_,__,C_,__,__,__,G_,__,__,__,__,__,__,N_,__, // 4x (upper case)
	__,__,__,__,T_,__,__,__,__,__,__,__,__,__,__,__, // 5x (upper case)
	__,a_,__,c_,__,__,__,g_,__,__,__,__,__,__,N_,__, // 6x (lower case)
	__,__,__,__,t_,__,__,__,__,__,__,__,__,__,__,__, // 7x (lower case)
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 8x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 9x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ax
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Bx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Cx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Dx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ex
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__	 // Fx
	};

static void do_colorize
   (u8*		dnaSeq,
	u8*		colorSeq,
	unspos	seqLen)
	{
	u8		bits1, bits2;

	bits1 = nuc_to_color_bits[X_];

	while (seqLen-- > 0)
		{
		bits2 = nuc_to_color_bits[*(dnaSeq++)];

		if (bits1 == Z_)
			*(colorSeq++) = 0;
		else if ((bits1 == X_) || (bits2 == X_))
			*(colorSeq++) = 'X';
		else if ((bits1 == N_) || (bits2 == N_))
			*(colorSeq++) = 'N';
		else if ((bits1 <= T_) || (bits2 <= T_))
			*(colorSeq++) = '0' + ((bits1 ^ bits2) & 3);
		else
			*(colorSeq++) = '4' + ((bits1 ^ bits2) & 3);

		bits1 = bits2;
		}

	*colorSeq = 0;
	}

#undef A_
#undef C_
#undef G_
#undef T_
#undef a_
#undef c_
#undef g_
#undef t_
#undef N_
#undef X_
#undef __
#undef Z_

//----------
//
// rev_comp_sequence--
//	Convert a sequence to its reverse-complement, in place.  Note that in
//	paratitioned sequences, each partition is reverse-complemented separately.
//
//----------
//
// Arguments:
//	seq*		_seq:				The sequence to reverse-complement.
//	const u8*	nucToComplement:	Array to map a nucleotide to its complement.
//									.. If this is NULL, nuc_to_complement is
//									.. used for DNA sequences.  For quantum
//									.. sequences this cannot be NULL.
//
// Returns:
//	(nothing)
//
//----------

static void revcomp_in_place (u8* seq, unspos seqLen, const u8* nucToComplement);

void rev_comp_sequence
   (seq*			_seq,
	const u8*		_nucToComplement)
	{
	seqpartition*	sp = &_seq->partition;
	partition*		p;
	u32				ix;
	const u8*		nucToComplement;

	if (_seq == NULL) suicide ("rev_comp_sequence(NULL)");
	if (_seq->len < 1) return;

	if ((_seq->fileType == seq_type_qdna)
	 && (_nucToComplement == NULL))
		{
		suicidef ("quantum DNA cannot be complemented (%s)\n"
		          "(the score file lacks complements)",
		          sequence_filename(_seq));
		return; // (we never reach here)
		}

	if (_nucToComplement != NULL) nucToComplement = _nucToComplement;
	                         else nucToComplement = nuc_to_complement;

	if (sp->p == NULL)
		revcomp_in_place (_seq->v, _seq->len, nucToComplement);
	else
		{
		p = sp->p;
		for (ix=0 ; ix<sp->len ; ix++)
			revcomp_in_place (/*start*/  _seq->v + p[ix].sepPos+1,
			                  /*length*/ p[ix+1].sepPos - (p[ix].sepPos+1),
			                  /*how*/    nucToComplement);
		}

	if (_seq->vc != NULL)
		do_colorize (_seq->v, _seq->vc, _seq->len);

	_seq->revCompFlags = _seq->revCompFlags ^ rcf_revcomp;
	}

static void revcomp_in_place
   (u8*			_seq,
	unspos		seqLen,
	const u8*	nucToComplement)
	{
	u8*		left, *right;
	u8		nuc;

	left  = _seq;
	right = left + seqLen-1;
	for ( ; left<=right ; left++,right--)
		{
		nuc    = nucToComplement[*left ];
		*left  = nucToComplement[*right];
		*right = nuc;
		}
	}

//----------
//
// backward_sequence--
//	Convert a sequence to its reverse (without complement), in place.  Note
//	that in paratitioned sequences, each partition is reversed separately.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to reverse.
//
// Returns:
//	(nothing)
//
//----------

static void backward_in_place (u8* seq, unspos seqLen);

void backward_sequence
   (seq*			_seq)
	{
	seqpartition*	sp = &_seq->partition;
	partition*		p;
	u32				ix;

	if (_seq->fileType == seq_type_csfasta)
		suicidef ("color space cannot be reversed (%s)", sequence_filename(_seq));

	if (_seq == NULL) suicide ("backward_sequence(NULL)");
	if (_seq->len < 1) return;

	if (sp->p == NULL)
		backward_in_place (_seq->v, _seq->len);
	else
		{
		p = sp->p;
		for (ix=0 ; ix<sp->len ; ix++)
			backward_in_place (/*start*/  _seq->v + p[ix].sepPos+1,
			                   /*length*/ p[ix+1].sepPos - (p[ix].sepPos+1));
		}

	_seq->revCompFlags = _seq->revCompFlags ^ rcf_rev;
	}


static void backward_in_place
   (u8*		_seq,
	unspos	seqLen)
	{
	u8*		left, *right;
	u8		nuc;

	left  = _seq;
	right = left + seqLen-1;
	for ( ; left<=right ; left++,right--)
		{ nuc = *left;  *left = *right;  *right = nuc; }
	}

//----------
//
// upper_sequence--
//	Convert a sequence to its upper-case equivalent, in place.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to modify.
//
// Returns:
//	(nothing)
//
//----------

static void upper_in_place (u8* seq, unspos seqLen);

void upper_sequence
   (seq*	_seq)
	{
	if (_seq == NULL) suicide ("upper_sequence(NULL)");
	if (_seq->len < 1) return;

	if (_seq->fileType == seq_type_csfasta)
		suicidef ("color space cannot be upper-cased (%s)", sequence_filename(_seq));
	else if (_seq->fileType == seq_type_qdna)
		suicidef ("quantum DNA cannot be upper-cased (%s)", sequence_filename(_seq));

	upper_in_place (_seq->v, _seq->len);
	}

static void upper_in_place
   (u8*		_seq,
	unspos	seqLen)
	{
	u8*		left, *right;

	left  = _seq;
	right = left + seqLen;
	for ( ; left<right ; left++)
		{
		if (('a' <= *left) && (*left <= 'z'))
			*left += 'A'-'a';
		}
	}

//----------
//
// parse_sequence_name--
//	Parse a sequence name
//
// The seqence name is the name of a file, plus some control options.  The
// basic form is
//
//	nickname::sequence_file/contig_name{mask_file}[actions]-
//
//----------
//
// Arguments:
//	const char*	name:				The name string to parse
//	char**		fileName:			Place to return the file name.
//	char**		nickname:			Place to return the nickname, if any.
//	char**		contigOfInterest:	Place to return the name of the particular
//									contig-of-interest, if any.
//	char**		namesFileName:		Place to return the contigs-of-interest
//									.. file name, if any.
//	int*		subsampleK:			Place to return K-of-N subsampling
//	int*		subsampleN:			.. specification.
//	char**		softMaskFileName:	Place to return the soft-mask file name, if
//									.. any.
//	int*		softMaskComplement:	Place to return true/false for soft-masking.
//	char**		xMaskFileName:		Place to return the x-mask file name, if any.
//	int*		xMaskComplement:	Place to return true/false for x-masking.
//	char**		nMaskFileName:		Place to return the n-mask file name, if any.
//	int*		xMaskComplement:	Place to return true/false for n-masking.
//	int*		nameParseType:		Place to return the name parse type, if any.
//	char**		nameTrigger:		Place to return the mask name trigger, if
//									.. any.
//	int*		doRevCompFlags:		Place to return rcf_forward/rcf_revcomp.
//									.. Actually can return any combination of
//									.. rcf_xxx flags.
//	int*		doUnmask, doJoin, useFullNames:
//									Place to return true/false for the
//									corresponding field in the sequence
//									structure.
// 	int*		fileType:			Place to return file type if any is
//									.. specified (one of seq_type_xxx).
//	int*		isQuantum:			Place to return true/false for whether the
//									.. sequence is to be quantum DNA.
//	char**		qCodingFileName:	Place to return the quantum coding file
//									.. name, if any.
//	unspos*		start:				Place to return the starting index (origin-1).
//									.. If no start is specified, 0 is placed
//									.. here. (see note 1 below)
//	unspos*		end:				Place to return the ending index (inclusive).
//									.. If no end is specified, 0 is placed here.
//									.. (see note 1 below)
//	int*		endIsSoft:			Place to return a flag telling the caller
//									.. that, while the end has been specified,
//									.. it is a soft specification and the caller
//									.. can trim it to the actual end of the
//									.. sequence.
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------
//
// Notes:
//
// (1)	Start and end are origin-1, inclusive on both ends.  So, for example,
//		start=10 and end=15 defines a 6 letter sequence, the 10th thru 15th
//		letters of the file.  The first 9 characters from the file should be
//		skipped.
//
//----------

void print_file_actions (FILE* f)
	{
	int   typeIx, indent, nameLen, needComma, lineWidth;
	char* action, *description, *name;

	fprintf (f, "Supported actions:\n");
	fprintf (f, "  <subrange>          only process a subrange of the file (see below)\n");
	fprintf (f, "  revcomp             reverse complement\n");
	fprintf (f, "  multiple            file's sequences are internally treated as a single\n");
	fprintf (f, "                      sequence\n");
	fprintf (f, "  subset=<namesfile>  process only the sequences listed in namesfile\n");
	fprintf (f, "                      (only valid for fasta, 2bit and hsx)\n");
	fprintf (f, "  subsample=<k>/<n>   process only the kth sequence of every group of n\n");
	fprintf (f, "                      sequences.  k ranges from 1 to n\n");
	fprintf (f, "                      (only valid for fasta, 2bit and hsx)\n");
	fprintf (f, "  unmask              convert any lowercase bases to uppercase\n");
	fprintf (f, "  softmask=<file>     mask segments specified in <file>, replacing them with\n");
	fprintf (f, "                      lowercase equivalents\n");
	fprintf (f, "  softmask=keep:<file> mask bases NOT in segments specified in <file>, with Xs\n");
	fprintf (f, "  xmask=<file>        mask segments specified in <file>, replacing them with Xs\n");
	fprintf (f, "  xmask=keep:<file>   mask bases NOT in segments specified in <file>, with Xs\n");
	fprintf (f, "  nmask=<file>        mask segments specified in <file>, replacing them with Ns\n");
	fprintf (f, "  nmask=keep:<file>   mask bases NOT in segments specified in <file>, with Ns\n");
	fprintf (f, "  nickname=<name>     name to use for this sequence in any output files\n");
	fprintf (f, "  nameparse=full      report full names in alignments instead of short names\n");
	fprintf (f, "  nameparse=alphanum  pull short name from sequence header, alphanumeric only\n");
	fprintf (f, "  nameparse=darkspace pull short name from sequence header, non-whitespace only\n");
	fprintf (f, "  nameparse=tag:<marker> pull a short name from sequence header, starting from\n");
	fprintf (f, "                      marker (only valid for fasta)\n");
	fprintf (f, "  quantum             the sequence contains quantum DNA\n");
	fprintf (f, "  quantum=<codefile>  the sequence contains quantum DNA, and <codefile>\n");
	fprintf (f, "                      describes the mapping from symbols to probabilities (only\n");
	fprintf (f, "                      meaningful for --format=text)\n");

	action      = "  format=<type>       ";
	description = "override auto-format detect;  <type> is one of ";
	fprintf (f, "%s%s", action, description);
	indent    = strlen (action);
	lineWidth = indent + strlen (description);
	for (typeIx=seq_type_unknown+1 ; typeIx<seq_type_max ; typeIx++)
		{
		name    = seqTypeNames[typeIx];
		nameLen = strlen(name);

		needComma = true;
		if (typeIx == seq_type_unknown+1)
			needComma = false;
		else if (lineWidth + 2 + nameLen >= 79)
			{
			fprintf (f, ",\n%*s", indent, " ");
			lineWidth = indent;
			needComma = false;
			}

		if (needComma) fprintf (f, ", ");
		fprintf (f, "%s", name);
		if (needComma) lineWidth += 2;
		lineWidth += nameLen;
		}
	fprintf (f, "\n\n");

	fprintf (f, "Subranges:\n");
	fprintf (f, "  start,end           process from start thru end, inclusive\n");
	fprintf (f, "  start..end          process from start thru end, inclusive\n");
	fprintf (f, "  start..             process from given start thru the end of the sequence\n");
	fprintf (f, "  ..end               process from the start of the sequence thru given end\n");
	fprintf (f, "  start#length        same as start..start+length-1\n");
	fprintf (f, "  center^length       same as center-length/2..center+length/2-1\n");
	fprintf (f, "  (subrange indices begin with 1 and are inclusive)\n");
	}


static void parse_sequence_name
   (const char*	name,
	char**		fileName,
	char**		nickname,
	char**		contigOfInterest,
	char**		namesFileName,
	int*		subsampleK,
	int*		subsampleN,
	char**		softMaskFileName,
	int*		softMaskComplement,
	char**		xMaskFileName,
	int*		xMaskComplement,
	char**		nMaskFileName,
	int*		nMaskComplement,
	int*		nameParseType,
	char**		nameTrigger,
	int*		doRevCompFlags,
	int*		doUnmask,
	int*		doJoin,
	int*		useFullNames,
	int*		fileType,
	int*		isQuantum,
	char**		qCodingFileName,
	unspos*		_start,
	unspos*		_end,
	int*		endIsSoft)
	{
	int			len;
	char*		fname, *bracket, *mask, *actions, *action, *actionName;
	char*		parse, *slashParse, *extParse;
	int			numItems, charsUsed;
	unspos		start, end, pendingStart, temp;
	int			tempInt;
	float		size;
	int			parsed;

	*namesFileName    = NULL;
	*contigOfInterest = NULL;
	*qCodingFileName  = NULL;

	*_start = *_end = 0;
	*endIsSoft = false;

	*doRevCompFlags = rcf_forward;
	*doUnmask       = false;
	*doJoin         = false;
	*useFullNames   = false;
	*fileType       = seq_type_unknown;
	*isQuantum      = false;

	//////////
	// copy the name, splitting out the nickname if present;  we will shorten
	// this copy if other components are present, so we are potentially
	// allocating more memory for the copy than is really needed
	//////////

	if (name == NULL) suicide ("parse_sequence_name(NULL)");

	parse   = strstr (name, "::");
	actions = strchr (name, '[');
	if ((parse == NULL)								// no "::"
	 || ((actions != NULL) && (parse > actions)))	// "::" is after "["
		{
		*nickname = NULL;
		*fileName = fname = copy_string (name);
		}
	else
		{
		if (parse-name == 0) goto empty_species_name;
		*nickname = copy_prefix (name, parse-name);
		*fileName = fname = copy_string (parse+2);
		}

	len = strlen (fname);
	if (len < 1) goto empty_file_name;

	//////////
	// see if we are to reverse the sequence
	//////////

	switch (fname[len-1])
		{
		case '-': *doRevCompFlags = rcf_revcomp;  fname[--len] = 0;  break;
		case '+':                                 fname[--len] = 0;  break;
		}

	if (len < 1) goto empty_file_name;

	//////////
	// split the file name string into its components
	//////////

	mask    = strchr (fname, '{');
	bracket = strchr (fname, '[');
	if ((bracket != NULL) && (mask != NULL))
		{ if (mask > bracket) mask = NULL; }

	if (mask == fname) goto empty_file_name;

	if (mask != NULL)
		*(mask++) = 0;

	parse = (mask == NULL)? fname : mask;

	actions = strchr (parse, '[');
	if (actions == parse) goto empty_file_name;

	if (actions != NULL)
		*(actions++) = 0;

	//////////
	// parse the mask file name
	//////////

	*softMaskFileName = NULL;
	*xMaskFileName    = NULL;
	*nMaskFileName    = NULL;

	if (mask != NULL)
		{
		len = strlen (mask);
		if (mask[--len] != '}') goto bad_mask;
		if (len == 0)           goto empty_mask_file_name;
		mask[len] = 0;
		*xMaskFileName = copy_string (mask);
		}

	//////////
	// split out the contig-of-interest if present
	//////////

	slashParse = NULL;

	extParse = strstr (fname, ".2bit/");
	if (extParse != NULL)
		slashParse = extParse+5;
	else
		{
		extParse = strstr (fname, ".hsx/");
		if (extParse != NULL)
			slashParse = extParse+4;
		}

	if (slashParse != NULL)
		{
		if      (strchr (slashParse+1, pathSlash) != NULL) extParse = NULL;
		else if (strstr (slashParse+1, ".2bit")   != NULL) extParse = NULL;
		}

	if (extParse != NULL)
		{
		*contigOfInterest = copy_string (slashParse+1);
		*slashParse = 0;
		}

	//////////
	// parse the actions list
	//////////

	if (actions != NULL)
		{
		len = strlen (actions);
		if (len == 0)
			goto bad_action_list;
		else if (actions[len-1] != ']')
			{
			if (strchr (actions, ']') == NULL) goto bad_action_list;
			                              else goto actions_not_at_end;
			}
		actions[--len] = 0;

		start = end = pendingStart = 0;
		while (actions != NULL)
			{
			//fprintf(stderr,"actions=\"%s\"\n", actions);
			action  = actions;
			actions = strchr (actions, ',');
			if (actions != NULL)
				*(actions++) = 0;
			else
				{
				actions = strstr (action, "][");
				if (actions != NULL)
					{ *(actions++) = 0;  actions++; }
				}

			//fprintf(stderr,"  action=\"%s\"\n", action);
			len = strlen (action);
			if (len == 0) goto blank_action;

			// parse simple actions

			if (strcmp (action, "unmask") == 0)
				{
				if (pendingStart != 0) goto unsatisfied_start;
				*doUnmask = true;
				continue;
				}

			if (strcmp (action, "revcomp") == 0)
				{
				if (pendingStart != 0) goto unsatisfied_start;
				*doRevCompFlags ^= rcf_revcomp;
				continue;
				}

			if ((strcmp (action, "backward")  == 0)  // (unadvertised)
			 || (strcmp (action, "backwards") == 0))
				{
				if (pendingStart != 0) goto unsatisfied_start;
				*doRevCompFlags ^= rcf_rev;
				continue;
				}

			if ((strcmp (action, "multi")    == 0)
			 || (strcmp (action, "multiple") == 0))
				{
				if (pendingStart != 0) goto unsatisfied_start;
				*doJoin = true;
				continue;
				}

			if ((strcmp (action, "nameparse=full") == 0)
			 || (strcmp (action, "fullname")       == 0)
			 || (strcmp (action, "fullnames")      == 0))
				{
				if (pendingStart != 0) goto unsatisfied_start;
				*useFullNames = true;
				continue;
				}

			if (action[0] == '@')
				{
				actionName = action + strlen("@");
				goto action_subset;
				}

			if (strcmp_prefix (action, "subset=") == 0)
				{
				actionName = action + strlen("subset=");
			action_subset:
				if (pendingStart != 0) goto unsatisfied_start;
				if (*namesFileName != NULL) goto many_name_files;
				if (strlen(actionName) == 0) goto bad_name_file;
				*namesFileName = copy_string (actionName);
				continue;
				}

			if (strcmp_prefix (action, "subsample=") == 0)
				{
				actionName = action + strlen("subsample=");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*subsampleN != 0) goto many_subsamples;
				if (strlen(actionName) == 0) goto bad_subsample;
				slashParse = strchr (actionName, '/');
				if (slashParse == NULL) goto bad_subsample;

				len = slashParse - actionName;
				*slashParse = ']';	// (write a sentinel for parsing K)
				charsUsed = -1;
				numItems = sscanf (actionName, "%d]%n", &tempInt, &charsUsed);
				if ((numItems != 1) || (charsUsed != len+1) || (tempInt < 1))
					{ *slashParse = '/';  goto bad_subsample; }
				*subsampleK = tempInt;

				*(slashParse++) = '/';
				len = strlen(slashParse);
				slashParse[len] = ']';	// (write a sentinel for parsing N)
				charsUsed = -1;
				numItems = sscanf (slashParse, "%d]%n", &tempInt, &charsUsed);
				if ((numItems != 1) || (charsUsed != len+1) || (tempInt < *subsampleK))
					{ slashParse[len] = 0;  goto bad_subsample; }
				*subsampleN = tempInt;
				continue;
				}

			if ((strcmp (action, "nameparse=alnum")    == 0)
			 || (strcmp (action, "nameparse=alphanum") == 0)
			 || (strcmp (action, "name:alnum")         == 0)
			 || (strcmp (action, "name:alphanum")      == 0))
				{
				if (pendingStart != 0) goto unsatisfied_start;
				if (*nameTrigger != NULL) goto many_name_parse_types;
				*nameParseType = name_parse_type_alnum;
				continue;
				}

			if (strcmp (action, "nameparse=darkspace") == 0)
				{
				if (pendingStart != 0) goto unsatisfied_start;
				if (*nameTrigger != NULL) goto many_name_parse_types;
				*nameParseType = name_parse_type_darkspace;
				continue;
				}

			if (strcmp_prefix (action, "nickname=") == 0)
				{
				actionName = action + strlen("nickname=");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*nickname != NULL)  goto many_nicknames;
				if (strlen(actionName) == 0) goto bad_nickname;
				*nickname = copy_string (actionName);
				continue;
				}

			if (strcmp_prefix (action, "name=") == 0)
				{
				actionName = action + strlen("name=");
				goto action_tag;
				}

			if (strcmp_prefix (action, "nameparse=tag:") == 0)
				{
				actionName = action + strlen("nameparse=tag:");
			action_tag:
				if (pendingStart != 0) goto unsatisfied_start;
				if (*nameTrigger != NULL) goto many_name_triggers;
				if (strlen(actionName) == 0) goto bad_name_trigger;
				if (*nameParseType != name_parse_type_core) goto many_name_parse_types;
				*nameParseType = name_parse_type_trigger;
				*nameTrigger   = copy_string (actionName);
				continue;
				}

			if (strcmp_prefix (action, "soft=keep:") == 0)
				{
				actionName = action + strlen("soft=keep:");
				goto action_softmaskkeep;
				}

			if (strcmp_prefix (action, "softmask=keep:") == 0)
				{
				actionName = action + strlen("softmask=keep:");
			action_softmaskkeep:
				if (pendingStart != 0) goto unsatisfied_start;
				if (*softMaskFileName != NULL) goto many_soft_mask_files;
				if (strlen(actionName) == 0) goto bad_soft_mask_file;
				*softMaskFileName   = copy_string (actionName);
				*softMaskComplement = true;
				continue;
				}

			if (strcmp_prefix (action, "soft=") == 0)
				{
				actionName = action + strlen("soft=");
				goto action_softmask;
				}

			if (strcmp_prefix (action, "softmask=") == 0)
				{
				actionName = action + strlen("softmask=");
			action_softmask:
				if (pendingStart != 0) goto unsatisfied_start;
				if (*softMaskFileName != NULL) goto many_soft_mask_files;
				if (strlen(actionName) == 0) goto bad_soft_mask_file;
				*softMaskFileName   = copy_string (actionName);
				*softMaskComplement = false;
				continue;
				}

			if (strcmp_prefix (action, "xmask=keep:") == 0)
				{
				actionName = action + strlen("xmask=keep:");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*xMaskFileName != NULL) goto many_x_mask_files;
				if (strlen(actionName) == 0) goto bad_x_mask_file;
				*xMaskFileName   = copy_string (actionName);
				*xMaskComplement = true;
				continue;
				}

			if (strcmp_prefix (action, "xmask=") == 0)
				{
				actionName = action + strlen("xmask=");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*xMaskFileName != NULL) goto many_x_mask_files;
				if (strlen(actionName) == 0) goto bad_x_mask_file;
				*xMaskFileName   = copy_string (actionName);
				*xMaskComplement = false;
				continue;
				}

			if (strcmp_prefix (action, "nmask=keep:") == 0)
				{
				actionName = action + strlen("nmask=keep:");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*nMaskFileName != NULL) goto many_n_mask_files;
				if (strlen(actionName) == 0) goto bad_n_mask_file;
				*nMaskFileName   = copy_string (actionName);
				*nMaskComplement = true;
				continue;
				}

			if (strcmp_prefix (action, "nmask=") == 0)
				{
				actionName = action + strlen("nmask=");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*nMaskFileName != NULL) goto many_n_mask_files;
				if (strlen(actionName) == 0) goto bad_n_mask_file;
				*nMaskFileName   = copy_string (actionName);
				*nMaskComplement = false;
				continue;
				}

			if (strcmp (action, "quantum") == 0)
				{
				if (pendingStart != 0) goto unsatisfied_start;
				if (*isQuantum) goto many_quantums;
				if (*fileType != seq_type_unknown) goto many_file_types;
				*isQuantum = true;
				*fileType = seq_type_qdna;
				continue;
				}

			if (strcmp_prefix (action, "quantum=") == 0)
				{
				actionName = action + strlen("quantum=");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*isQuantum) goto many_quantums;
				if (strlen(actionName) == 0) goto bad_code_file;
				*qCodingFileName = copy_string (actionName);
				*isQuantum       = true;
				continue;
				}

			if (strcmp_prefix (action, "format=") == 0)
				{
				int fType, typeIx;
				actionName = action + strlen("format=");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*fileType != seq_type_unknown) goto many_file_types;
				fType = seq_type_unknown;
				for (typeIx=seq_type_unknown+1 ; typeIx<seq_type_max ; typeIx++)
					{
					if (strcmp (actionName, seqTypeNames[typeIx]) == 0)
						{ fType = typeIx;  break; }
					}
				if (fType == seq_type_unknown) goto bad_file_format;
				*fileType = typeIx;
				continue;
				}

			// any other option would be a sequence limit;  so if we already
			// have sequence limits, this is an error

			if ((start != 0) || (end != 0))
				goto bad_action;

			// try to parse the field as a single number

			action[len] = ']';		// (write a sentinel;  this location is
									//  guaranteed to be part of the string)

			charsUsed = -1;
			numItems = sscanf (action, unsposFmtScanf "]%n", &temp, &charsUsed);
			if ((numItems == 1) && (charsUsed == len+1))
				{
				//fprintf(stderr,"  parsed as one (" unsposFmt ") %d\n", temp, len);
				if (temp == 0)
					{ action[len] = 0;  goto bad_sequence_position; }
				if (pendingStart == 0)
					pendingStart = temp;
				else
					{ start = pendingStart;  end = temp;  pendingStart = 0; }
				continue;
				}
			else if (pendingStart != 0)
				{ action[len] = 0;  goto unsatisfied_start; }

			// try to parse the field as sequence limits

			parsed = false;

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposDotsFmtScanf "]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as two, dots\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposCommaFmtScanf "]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as two, comma\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "#" unsposFmtScanf "]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as two, waffle\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end += start-1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "#" unsposFmtScanf "K]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, waffle with K\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end *= 1000;
					end += start-1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "#%fK]%n", &start, &size, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, waffle with K\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end =  (size * 1000) + 1;
					end += start-1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "#" unsposFmtScanf "M]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, waffle with M\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end *= 1000 * 1000;
					end += start-1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "#%fM]%n", &start, &size, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, waffle with M\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end =  (size * 1000 * 1000) + 1;
					end += start-1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "^" unsposFmtScanf "]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as two, caret\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					start -= end / 2;
					end   += start-1;
					if (start < 1) start = 1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "^" unsposFmtScanf "K]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, caret with K\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end   *= 1000;
					start -= end / 2;
					end   += start-1;
					if (start < 1) start = 1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "^%fK]%n", &start, &size, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, caret with K\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end   =  (size * 1000) + 1;
					start -= end / 2;
					end   += start-1;
					if (start < 1) start = 1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "^" unsposFmtScanf "M]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, caret with M\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end   *= 1000 * 1000;
					start -= end / 2;
					end   += start-1;
					if (start < 1) start = 1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "^%fM]%n", &start, &size, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, caret with M\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end   =  (size * 1000 * 1000) + 1;
					start -= end / 2;
					end   += start-1;
					if (start < 1) start = 1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "..]%n", &start, &charsUsed);
				if ((numItems == 1) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, dots at end\n");
					if (start == 0)
						{ action[len] = 0;  goto bad_limits; }
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, ".." unsposFmtScanf "]%n", &end, &charsUsed);
				if ((numItems == 1) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, dots at start\n");
					if (end == 0)
						{ action[len] = 0;  goto bad_limits; }
					parsed = true;
					}
				}

			if (!parsed)
				{
				action[len] = 0;		// (clear the sentinel)
				goto bad_action;
				}
			}

		if (pendingStart != 0) goto unsatisfied_start_2;

		if ((start != 0) && (end != 0) && (start > end))
			{
			temp = start;  start = end;  end = temp;
			*doRevCompFlags ^= rcf_revcomp;
			}

		if ((start != 0) || (end != 0))
			{
			*_start = start;
			*_end   = end;
			}
		}

	//////////
	// sanity checks
	//////////

	if ((*contigOfInterest != NULL) && (*namesFileName != NULL))
		suicidef ("(for %s) can't use these together:"
		          "         %s\n"
		          "         %s",
		          fname, contigOfInterest, namesFileName);

	if ((*doJoin) && (*softMaskFileName != NULL))
		suicidef ("(for %s) can't use multi with softmask=%s",
		          fname, softMaskFileName);

	if ((*doJoin) && (*xMaskFileName != NULL))
		suicidef ("(for %s) can't use multi with xmask=%s",
		          fname, xMaskFileName);

	if ((*doJoin) && (*nMaskFileName != NULL))
		suicidef ("(for %s) can't use multi with nmask=%s",
		          fname, nMaskFileName);

	return;

	//////////
	// failure exits
	//////////

empty_file_name:
	suicidef ("sequence file name is absent from \"%s\"", name);
	return; // (can't reach here)

empty_species_name:
	suicidef ("(for %s) empty nickname", parse+2);
	return; // (can't reach here)

bad_mask:
	suicidef ("(for %s) mask name needs closing } (%s)", fname, mask);
	return; // (can't reach here)

empty_mask_file_name:
	suicidef ("(for %s) use [unmask] instead of {}", fname);
	return; // (can't reach here)

actions_not_at_end:
	suicidef ("(for %s) action list is not at end", fname);
	return; // (can't reach here)

bad_action_list:
	suicidef ("(for %s) bad action list", fname);
	return; // (can't reach here)

blank_action:
	suicidef ("(for %s) blank action", fname);
	return; // (can't reach here)

bad_action:
	suicidef ("(for %s) bad action \"%s\"", fname, action);
	return; // (can't reach here)

bad_sequence_position:
	suicidef ("(for %s) bad limit \"%s\"", fname, action);
	return; // (can't reach here)

bad_limits:
	suicidef ("(for %s) bad limits \"%s\"", fname, action);
	return; // (can't reach here)

unsatisfied_start:
	suicidef ("(for %s) incomplete limits (%d,%s)", fname, pendingStart, action);
	return; // (can't reach here)

unsatisfied_start_2:
	suicidef ("(for %s) incomplete limits (%d)", fname, pendingStart);
	return; // (can't reach here)

many_name_files:
	suicidef ("(for %s) only one names file allowed:\n"
	          "         subset=%s\n"
	          "         subset=%s",
	          fname, *namesFileName, actionName);
	return; // (can't reach here)

bad_name_file:
	suicidef ("(for %s) subset= requires a names file", fname);
	return; // (can't reach here)

many_subsamples:
	suicidef ("(for %s) only one subsampling allowed:\n"
	          "         subsample=%d/%d\n"
	          "         subsample=%s",
	          fname, *subsampleK, *subsampleN, actionName);
	return; // (can't reach here)

bad_subsample:
	suicidef ("(for %s) bad subsample \"%s\"", fname, actionName);
	return; // (can't reach here)

many_nicknames:
	suicidef ("(for %s) only one nickname allowed:\n"
	          "         nickname=%s\n"
	          "         nickname=%s",
	          fname, *nickname, actionName);
	return; // (can't reach here)

bad_nickname:
	suicidef ("(for %s) nickname= requires a non-empty string", fname);
	return; // (can't reach here)

many_name_parse_types:
	suicidef ("(for %s) only one name parsing allowed\n", fname);
	return; // (can't reach here)

many_name_triggers:
	suicidef ("(for %s) only one name trigger allowed:\n"
	          "         nameparse=tag:%s\n"
	          "         nameparse=tag:%s",
	          fname, *nameTrigger, actionName);
	return; // (can't reach here)

bad_name_trigger:
	suicidef ("(for %s) nameparse=tag: requires a non-empty string", fname);
	return; // (can't reach here)

many_soft_mask_files:
	suicidef ("(for %s) only one softmask allowed:\n"
	          "         softmask=%s\n"
	          "         softmask=%s",
	          fname, *softMaskFileName, actionName);
	return; // (can't reach here)

bad_soft_mask_file:
	suicidef ("(for %s) softMask= or softMask=keep: require a non-empty string", fname);
	return; // (can't reach here)

many_x_mask_files:
	suicidef ("(for %s) only one xmask allowed:\n"
	          "         xmask=%s\n"
	          "         xmask=%s",
	          fname, *xMaskFileName, actionName);
	return; // (can't reach here)

bad_x_mask_file:
	suicidef ("(for %s) xmask= or xmask=keep: require a non-empty string", fname);
	return; // (can't reach here)

many_n_mask_files:
	suicidef ("(for %s) only one nmask allowed:\n"
	          "         nmask=%s\n"
	          "         nmask=%s",
	          fname, *nMaskFileName, actionName);
	return; // (can't reach here)

bad_n_mask_file:
	suicidef ("(for %s) nmask= or nmask=keep: require a non-empty string", fname);
	return; // (can't reach here)

many_file_types:
	suicidef ("(for %s) more than one file type is defined", fname);
	return; // (can't reach here)

bad_file_format:
	suicidef ("(for %s) unknown file format: %s", fname, actionName);
	return; // (can't reach here)

many_quantums:
	suicidef ("(for %s) only one instance of quantum allowed", fname);
	return; // (can't reach here)

bad_code_file:
	suicidef ("(for %s) quantum= requires a non-empty string", fname);
	return; // (can't reach here)
	}

//----------
//
// detect_file_type--
//	Attempt to determine what type of file we are dealing with (e.g. fasta,
//	nib, quantum, etc.).
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence.
//
// Returns:
//	The type of file being read (one of seq_type_xxx);  failure causes program
//	fatality.
//
//----------

static int detect_file_type
   (seq*	_seq)
	{
	int		type = seq_type_unknown;
	char	buffer[maxSequenceHeader+3];
	u32		bufferLen;
	u32		magic;
	u8		ch;
	int		intCh;

	//////////
	// determine if it's a nib file (from the magic number)
	//////////

	// read the first four bytes;  if it's a recognizable magic number then it
	// must be the corresponding type of file;  otherwise we assume it must be
	// a fasta file (if not, it will die later)

	magic = read_4_little (_seq);

	if ((magic == nibMagicLittle)          || (magic == nibMagicBig))
		type = seq_type_nib;
	else if ((magic == twobitMagicLittle)  || (magic == twobitMagicBig))
		type = seq_type_2bit;
	else if ((magic == hsxMagicLittle)     || (magic == hsxMagicBig))
		type = seq_type_hsx;
	else if ((magic == qdnaMagicLittle)    || (magic == qdnaMagicBig))
		type = seq_type_qdna;
	else if ((magic == oldQdnaMagicLittle) || (magic == oldQdnaMagicBig))
		type = seq_type_qdna;

	// put those four bytes back in the file (in reverse of the read order)

	seq_ungetc (magic >> 24, _seq);
	seq_ungetc (magic >> 16, _seq);
	seq_ungetc (magic >>  8, _seq);
	seq_ungetc (magic      , _seq);

	if (type != seq_type_unknown)
		return type;

	//////////
	// determine if it's a fasta or csfasta file;  these are very similar
	// formats;  if the first character is a '#' we know it is csfasta (because
	// regular fasta doesn't allow comments);  if the first character is a '>'
	// then it can still be either, in which case we need to skip the header
	// line and read the first two sequence characters
	//////////

	ch = seq_getc (_seq);
	if (ch == '#')
		{
		seq_ungetc (ch, _seq);
		return seq_type_csfasta;
		}

	if (ch != '>')
		seq_ungetc (ch, _seq);
	else
		{
		// read header

		bufferLen = 0;
		buffer[bufferLen++] = ch;

		do
			{
			intCh = seq_getc (_seq);
			if (intCh == EOF) goto unknown;
			if (bufferLen >= sizeof(buffer)-3) goto unknown;
			buffer[bufferLen++] = (char) intCh;
			} while (intCh != '\n');

		// first character must be a nucleotide or it is not fasta
		// $$$ we are ignoring the slim chance of an empty first sequence, or
		// $$$ .. that the first line contains just a single nucleotide

		intCh = seq_getc (_seq);
		if (intCh == EOF) goto unknown;
		buffer[bufferLen++] = intCh;
		if (ustrchr ("ACGTacgtNn", intCh) != NULL)
			{
			intCh = seq_getc (_seq);
			if (intCh == EOF) goto unknown;
			buffer[bufferLen++] = intCh;
			if (ustrchr ("ACGTacgtNn", intCh) != NULL)
				type = seq_type_fasta;
			else if (ustrchr ("0123", intCh) != NULL)
				type = seq_type_csfasta;
			}

	unknown:
		while (bufferLen > 0)
			seq_ungetc (buffer[--bufferLen], _seq);

		if (type != seq_type_unknown)
			return type;
		}

	//////////
	// if all else fails, assume it's a fasta file (for compatibility with
	// blastz)
	//////////

	if (type == seq_type_unknown)
		type = seq_type_fasta;

	return type;
	}

//----------
//
// read_4, read_4_big, read_4_little--
//	Read four bytes from a file, in big or little endian order.
//
//----------
//
// Arguments:
//	seq*	_seq:			The sequence.
//	int		asBigEndian:	true  => read 'em as big endian
//							false => read 'em as little endian
//
// Returns:
//	The magic number read.
//
//----------

static u32 read_4
   (seq*	_seq,
	int		asBigEndian)
	{
	if (asBigEndian) return read_4_big    (_seq);
	            else return read_4_little (_seq);
	}

static u32 read_4_big
   (seq*	_seq)
	{
	u32		val;

	val =  seq_getc (_seq) << 24;
	val |= seq_getc (_seq) << 16;
	val |= seq_getc (_seq) << 8;
	val |= seq_getc (_seq);

	return val;
	}

static u32 read_4_little
   (seq*	_seq)
	{
	u32		val;

	val =  seq_getc (_seq);
	val |= seq_getc (_seq) << 8;
	val |= seq_getc (_seq) << 16;
	val |= seq_getc (_seq) << 24;

	return val;
	}


static u64 read_5
   (seq*	_seq,
	int		asBigEndian)
	{
	if (asBigEndian) return read_5_big    (_seq);
	            else return read_5_little (_seq);
	}

static u64 read_5_big
   (seq*	_seq)
	{
	u64		val;

	val =  ((u64) seq_getc (_seq)) << 32;
	val |=        seq_getc (_seq)  << 24;
	val |=        seq_getc (_seq)  << 16;
	val |=        seq_getc (_seq)  << 8;
	val |=        seq_getc (_seq);

	return val;
	}

static u64 read_5_little
   (seq*	_seq)
	{
	u64		val;

	val =         seq_getc (_seq);
	val |=        seq_getc (_seq)  << 8;
	val |=        seq_getc (_seq)  << 16;
	val |=        seq_getc (_seq)  << 24;
	val |= ((u64) seq_getc (_seq)) << 32;

	return val;
	}


static u64 read_6
   (seq*	_seq,
	int		asBigEndian)
	{
	if (asBigEndian) return read_6_big    (_seq);
	            else return read_6_little (_seq);
	}

static u64 read_6_big
   (seq*	_seq)
	{
	u64		val;

	val =  ((u64) seq_getc (_seq)) << 40;
	val |= ((u64) seq_getc (_seq)) << 32;
	val |=        seq_getc (_seq)  << 24;
	val |=        seq_getc (_seq)  << 16;
	val |=        seq_getc (_seq)  << 8;
	val |=        seq_getc (_seq);

	return val;
	}

static u64 read_6_little
   (seq*	_seq)
	{
	u64		val;

	val =         seq_getc (_seq);
	val |=        seq_getc (_seq)  << 8;
	val |=        seq_getc (_seq)  << 16;
	val |=        seq_getc (_seq)  << 24;
	val |= ((u64) seq_getc (_seq)) << 32;
	val |= ((u64) seq_getc (_seq)) << 40;

	return val;
	}

//----------
//
// skip_seq_whitespace, skip_fasta_whitespace--
//	Read characters from the associated sequence until we get something that
//	ain't whitespace.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to read.
//
// Returns:
//	(same as for getc())
//
//----------

static int skip_seq_whitespace
   (seq*	_seq)
	{
	int		ch;

	do
		{
		ch = seq_getc (_seq);
		} while ((ch == ' ') || (ch == '\t'));

	return ch;
	}

static int skip_fasta_whitespace
   (seq*	_seq)
	{
	int		ch;

	do
		{
		ch = seq_getc (_seq);
		} while ((ch == ' ') || (ch == '\t'));

	return ch;
	}

//----------
//
// seq_getc--
//	Read the next character from the associated file.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to read.
//
// Returns:
//	(same as for getc();  the character read, or EOF)
//
//----------

static int seq_getc
   (seq*	_seq)
	{
	// if there are characters pending, get one from the pending buffer;
	// otherwise, feed one straight from the file

	if (_seq->pendingLen == 0)
		return getc_or_die (_seq->f, _seq->fileName);

	_seq->pendingLen--;
	return (int) (u8) *(_seq->pendingStack++);
	}

//----------
//
// seq_ungetc--
//	Give back a character to the associated file.
//
// WARNING: This routine is not an exact drop in replacement for the standard
//	c routine ungetc().
//
//----------
//
// Arguments:
//	char	ch:		The character to return.  Characters should be returned in
//					.. the opposite order of that read (ie newest char returned
//					.. first).
//	seq*	_seq:	The sequence being read.
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static void seq_ungetc
   (char	ch,
	seq*	_seq)
	{
	if (_seq->pendingLen >= seqBufferSize)
		suicide ("seq_ungetc() buffer is already full");

	_seq->pendingLen++;
	*(--_seq->pendingStack) = ch;
	}

//----------
//
// skip_chars--
//	Skip the next so many characters from the associated file.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence being read.
//	u32		toSkip:	The number of characters to skip.
//
// Returns:
//	true if successful, false if there's a problem (such as reaching premature
//	end-of-file).
//
//----------

static int skip_chars
   (seq*	_seq,
	u32		toSkip)
	{
	int		ch;

	// see if we have any to skip in the pending buffer

	if (_seq->pendingLen >= toSkip)
		{
		// we have all we need in the pending buffer

		_seq->pendingLen   -= toSkip;
		_seq->pendingStack += toSkip;
		return true;
		}

	if (_seq->pendingLen > 0)
		{
		// everything in the pending buffer will be skipped

		toSkip            -= _seq->pendingLen;
		_seq->pendingLen   =  0;
		_seq->pendingStack = _seq->pendingChars + seqBufferSize;
		}

	if (toSkip == 0) return true; // (none left to skip)

	// skip the rest by seeking past the characters

	if (fseek (_seq->f, toSkip, SEEK_CUR) != 0)
		{
		// seek failed, so let's try reading instead

		while (toSkip-- > 0)
			{
			ch = getc (_seq->f);
			if ((ch == EOF) || (ferror (_seq->f)))
				return false;
			}
		}

	return true;
	}

//----------
//
// test_rewindability--
//	Test wheteher a sequence's underflying file is rewinable.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence being read.
//
// Returns:
//	An fseek error code;  zero indicates the underlying file is rewindable;
//	any other value indicates that it is not.
//
//----------

static int test_rewindability
   (seq*		_seq)
	{
	long int	savedFilePos;

	savedFilePos = ftell (_seq->f);
	return fseek (_seq->f, savedFilePos, SEEK_SET);
	}

//----------
//
// save_fstate, restore_fstate--
//	Save and restore the state of the associated file.  This lets the caller
//	read ahead in a file, then return to the original point.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence being read.
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static void save_fstate
   (seq*	_seq)
	{
	if (_seq == NULL) suicide ("save_fstate(NULL)");

	// save read head

	_seq->savedFilePos  = ftell (_seq->f) - _seq->pendingLen;
	_seq->hasSavedState = true;
	}

static void restore_fstate
   (seq*	_seq)
	{
	int		err;

	if (_seq == NULL)         suicide ("restore_fstate(NULL)");
	if (!_seq->hasSavedState) suicide ("restore_fstate(), no state saved");

	// restore read head

	err = fseek (_seq->f, _seq->savedFilePos, SEEK_SET);
	if (err != 0)
		suicidef_with_perror ("restore_fstate(), fseek returned %d", err);

	// restore pending character state

	_seq->pendingLen   = 0;
	_seq->pendingStack = _seq->pendingChars + seqBufferSize;
	}

//----------
//
// match_composition--
//	Count the number of matched DNA letter pairs in a gap free alignment.
//
//----------
//
// Arguments:
//	seq*	seq1:		The first sequence.
//	unspos 	pos1:		The subsequence start position in seq1 (origin-0).
//	seq*	seq2:		The second sequence.
//	unspos 	pos2:		The subsequence start position in seq2 (origin-0).
//	unspos 	length:		The length of the subsequence.
//	unspos	count[4][4]:Place to return the counts of each matched DNA letter
//						.. pair.  Indexing is as per nuc_to_bits.
//
// Returns:
//	nothing;  composition is returned in the count[][] array
//
//----------
//
// Note:  Masked (lowercase) bp do not contribute to the results.
//
//----------

void match_composition
   (seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length,
	unspos	count[4][4])
	{
	u8*		s1 = seq1->v + pos1;
	u8*		s2 = seq2->v + pos2;
	unspos	ix;
	int		r, c;

	for (r=0 ; r<4 ; r++)
		for (c=0 ; c<4 ; c++)
			count[r][c] = 0;

	for (ix=0 ; ix<length ; ix++)
		{
		r = upper_nuc_to_bits[*(s1++)];
		c = upper_nuc_to_bits[*(s2++)];
		if ((r >= 0) && (c >= 0)) 
			count[r][c]++;
		}
	}

//----------
//
// percent_identical--
//	Determine the percentage of bases that match in two subsequences.
//
//----------
//
// Arguments:
//	seq*	seq1:	The first sequence.
//	unspos 	pos1:	The subsequence start position in seq1 (origin-0).
//	seq*	seq2:	The second sequence.
//	unspos 	pos2:	The subsequence start position in seq2 (origin-0).
//	unspos 	length:	The length of the subsequence.
//
// Returns:
//	The percentage (an integer in the range 0..100).
//
//----------
//
// Note:  Masked (lowercase) bp *do* contribute to the results, but illegal
//        values like N do not (and they are not counted in the denominator
//        either).
//
//----------

int percent_identical
   (seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length)
	{
	u8*		s1 = seq1->v + pos1;
	u8*		s2 = seq2->v + pos2;
	s8		c1, c2;
	unspos	numMatches = 0;
	unspos	denom = 0;
	unspos	ix;

	if (length == 0)
		return 0;

	if ((seq1->fileType == seq_type_qdna)
	 || (seq2->fileType == seq_type_qdna))
		return 0;

	for (ix=0 ; ix<length ; ix++)
		{
		c1 = nuc_to_bits[*(s1++)];
		c2 = nuc_to_bits[*(s2++)];
		if ((c1 >= 0) && (c2 >= 0))
			{
			if (c1 == c2) numMatches++;
			denom++;
			}
		}

	if (denom == 0)
		return 0;
	else
		return (200*numMatches + denom) / (2*denom); // 100*numMatches/denom, rounded
	}

//----------
//
// score_match--
//	Determine the substitution score of aligned bases in two subsequences.
//
//----------
//
// Arguments:
//	scoreset*	scoring:	The scoring scheme to use.
//	seq*		seq1:		The first sequence.
//	unspos 		pos1:		The subsequence start position in seq1 (origin-0).
//	seq*		seq2:		The second sequence.
//	unspos 		pos2:		The subsequence start position in seq2 (origin-0).
//	unspos 		length:		The length of the subsequence.
//
// Returns:
//	The substitution score of the two subsequences.
//
//----------

score score_match
   (scoreset*	scoring,
	seq*		seq1,
	unspos		pos1,
	seq*		seq2,
	unspos		pos2,
	unspos		length)
	{
	u8*			s1 = seq1->v + pos1;
	u8*			s2 = seq2->v + pos2;
	u8*			stop = s1 + length;
	score		similarity = 0;

	if (length == 0)
		return 0;

	while (s1 < stop)
		similarity += scoring->sub[*s1++][*s2++];

	return similarity;
	}

//----------
//
// dump_aligned_nucleotides--
//	Dump the nucleotides (from each sequence) for a gap-free alignment.
//
//----------
//
// Arguments:
//	FILE*	f:		The file to print to.
//	seq*	seq1:	One sequence.
//	unspos	pos1:	The first aligned position in sequence 1.
//	seq*	seq2:	The other sequence.
//	unspos	pos2:	The first aligned position in sequence 2.
//	unspos	length:	The length of the alignment.
//
// Returns:
//	(nothing)
//
//----------

void dump_aligned_nucleotides
   (FILE*	f,
	seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length)
	{
	int		isRev1 = ((seq1->revCompFlags & rcf_rev) != 0);
	int		isRev2 = ((seq2->revCompFlags & rcf_rev) != 0);
	char*	start1 = (char*) seq1->v + pos1;
	char*	start2 = (char*) seq2->v + pos2;
	int		digits = 10;

	fprintf      (f, unsposStarFmt "%c:", digits, pos1+1, (isRev1)?'-':'+');
	print_prefix (f, (char*) seq1->v + pos1, (int) length);
	fprintf      (f, "\n");

	fprintf      (f, "%*s  ", digits, "");
	print_dna_similarities
	             (f, start1, start2, (int) length);
	fprintf      (f, "\n");

	fprintf      (f, unsposStarFmt "%c:", digits, pos2+1, (isRev2)?'-':'+');
	print_prefix (f, (char*) seq2->v + pos2, (int) length);
	fprintf      (f, "\n");
	}

//----------
//
// dump_sequence--
//	Write a sequence to a file (for debugging).
//
//----------
//
// Arguments:
//	FILE*	f:		The file to print to.
//	seq*	_seq:	The sequence to print.
//
// Returns:
//	(nothing)
//
//----------

static void dump_sequence
   (FILE*	f,
	seq*	_seq)
	{
	char	buffer[101];
	unspos	ix, start = 0;
	int		width;
	char	ch;
	int		bx;

	start = _seq->len;
	width = 1;
	while (start > 9) { start/=10;  width++; }

	bx = 0;
	for (ix=0 ; ix<_seq->len ; ix++)
		{
		ch = (char) _seq->v[ix];

		if (ch == 0)
			{
			if (bx > 0)
				{
				buffer[bx] = 0;
				fprintf (f, unsposStarFmt ": %s\n", width, start, buffer);
				}
			continue;
			}

		if (bx == sizeof(buffer)-1)
			{
			buffer[bx] = 0;
			fprintf (f, unsposStarFmt ": %s\n", width, start, buffer);
			bx = 0;
			}

		if (bx == 0) start = ix;
		buffer[bx++] = ch;
		}

	if (bx > 0)
		{
		buffer[bx] = 0;
		fprintf (f, unsposStarFmt ": %s\n", width, start, buffer);
		}
	}

