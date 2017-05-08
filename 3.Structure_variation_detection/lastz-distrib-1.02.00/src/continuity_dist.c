//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: continuity_dist.c
//
//----------
//
// continuity_dist--
//	Support for collecting the query continuity (or gap rate) distribution from
//	alignments.
//
// Notes:
//	continuity = 1/(1+gaprate)
//  gaprate    = (1/continuity)-1
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
#include <ctype.h>				// standard C upper/lower stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

#define  continuity_dist_owner	// (make this the owner of its globals)
#include "continuity_dist.h"	// interface to this module

//----------
//
// filter_aligns_by_continuity--
//	Filter a list of alignments, removing any alignment that has continuity
//	percentage outside of a specified range.
//
//----------
//
// Arguments:
//	alignel*	alignList:		The list of alignments to operate upon.
//	float		minContinuity,	The range of query continuity that we will
//				maxContinuity:	.. *keep*. These are values between 0 and 1.
//
// Returns:
//	A pointer to the list of remaining alignments.
//
//----------
//
// Notes:
//	(1)	The numerator is the number of alignment columns *not* containing gaps
//      and  the denominator is the number of total alignment columns.
//	(2)	Memory for alignments that don't make the cut is deallocated here.
//	(3)	The returned list of alignments is in the same order as the incoming
//		list.
//
//----------

alignel* filter_aligns_by_continuity
   (alignel*	alignList,
	float		minContinuity,
	float		maxContinuity)
	{
	alignel*	a, *next;
	alignel*	head, *prev;
	unspos		numer, denom;
	float		minCont, maxCont;

	// process each alignment, collecting a list of those that are long enough

	head = prev = NULL;
	for (a=alignList ; a!=NULL ; a=next)
		{
		next  = a->next;

		// if the alignment is too gappy, skip it

		alignment_continuity (a, &numer, &denom);

		minCont = denom * minContinuity;
		maxCont = denom * maxContinuity;

		if ((numer < minCont) || (numer > maxCont))
			{ // (unwanted alignment, discard it)
			free_if_valid ("filter_aligns_by_continuity a->script", a->script);
			free_if_valid ("filter_aligns_by_continuity a",         a);
			continue;
			}


		// this alignment is ok, add it to the end of the new list we're
		// building

		if (head == NULL) head = prev = a;
		             else { prev->next = a;  prev = a; }

		a->next = NULL;
		}

	return head;
	}

//----------
//
// alignment_continuity--
//	Compute the continuity of an gapped alignment block.  This is the number of
//	bases *not* aligned to gaps divided by the number of alignment columns.
//
//----------
//
// Arguments:
//	alignel*	a:				The alignment of interest.
//	unspos*		numer, denom:	Place to return the continuity fraction.  Note
//								.. that the returned denominator might be zero,
//								.. and that the fraction may be greater than 1.
//
// Returns:
//	(nothing)
//
//----------

void alignment_continuity
   (alignel*	a,
	unspos*		_numer,
	unspos*		_denom)
	{
	unspos		gapColumns, nonGapColumns;

	alignment_gap_rate (a, &gapColumns, &nonGapColumns);

	*_numer = nonGapColumns;
	*_denom = nonGapColumns + gapColumns;
	}

//----------
//
// alignment_gap_rate--
//	Compute the gap rate of an gapped alignment block.  This is the number of
//	bases aligned to gaps divided by the number of aligned bases.
//
//----------
//
// Arguments:
//	alignel*	a:				The alignment of interest.
//	unspos*		numer, denom:	Place to return the gap rate fraction.  Note
//								.. that the returned denominator might be zero,
//								.. and that the fraction may be greater than 1.
//
// Returns:
//	(nothing)
//
//----------

void alignment_gap_rate
   (alignel*	a,
	unspos*		_numer,
	unspos*		_denom)
	{
	unspos		beg1 = a->beg1;
	unspos		beg2 = a->beg2;
	unspos		height, width, i, j, prevI, prevJ;
	u32			opIx;
    unspos		run;
	unspos		denom, gappedBases;

	height = a->end1 - beg1 + 1;
	width  = a->end2 - beg2 + 1;

	denom = 0;
	opIx  = 0;
	for (i=j=0 ; (i< height)||(j<width) ; )
		{
		prevI = i;  prevJ = j;
		run = edit_script_run_of_subs (a->script, &opIx);
		i += run; j += run;

		denom += run;

		if ((i < height) || (j < width))
			edit_script_indel_len (a->script, &opIx, &i, &j);
		}

	if (denom == 0)
		{ *_numer = *_denom = 0;  return; }

	gappedBases = (height - denom) + (width - denom);

	*_numer = gappedBases;
	*_denom = denom;
	}

