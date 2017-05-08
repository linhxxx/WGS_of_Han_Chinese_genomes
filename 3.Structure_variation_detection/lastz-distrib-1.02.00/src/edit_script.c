//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: edit_script.c
//
//----------
//
// edit_script--
//	Support for representing alignments as a series of substitute, insert,
//	and delete.
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
#include "utilities.h"			// utility stuff

#define  edit_script_owner		// (make this the owner of its globals)
#include "edit_script.h"		// interface to this module

//----------
//
// data structures and types
//
//----------

// edit operations (type editop)
//
// each consists of a 2-bit code (for insert, delete or substitute) and a
// 30-bit repeat count

#define edit_op_operation(op)     ((op) & 0x3)
#define edit_op_repeat(op)        ((op) >> 2)
#define edit_op(op,rpt)           (((op) & 0x3) | ((rpt) << 2))
#define edit_op_add_repeat(eop,n) ((eop) + ((n) << 2))

#define maxEditopRepeat ((((u32)1)<<30)-1)

//----------
//
// prototypes for private functions
//
//----------

static void edit_script_make_room (editscript** s, u32 entries);
static void edit_script_put       (editscript** s, u32 op, u32 rpt);

//----------
//
// free_align_list--
//	Dispose of a list of alignments.
//
//----------
//
// Arguments:
//	alignel*	a:	The list of alignments to dispose of.
//
// Returns:
//	(nothing)
//
//----------

void free_align_list
   (alignel*	a)
	{
	alignel*	b;

	while (a != NULL)
		{
		b = a->next;
		free_if_valid ("free_align_list a->script", a->script);
		free_if_valid ("free_align_list a",         a);
		a = b;
		}
	}

//----------
//
// edit_script_new--
//	Allocate an alignment edit script.
//
//----------
//
// Arguments:
//	(none)
//
// Returns:
//	A pointer to an empty edit script.  The caller is responsible for disposing
//	of this memory, for which purpose free() can be used.
//
//----------
//
// Notes:  NULL is never returned-- failure to allocate is a fatal error.
//
//----------

editscript* edit_script_new
   (void)
	{
	u32			entries = 12;
	editscript*	s;

	// allocate

	s = zalloc_or_die ("edit_script_new", edit_script_bytes(entries));

	// initialize;  note that by use of zalloc we already have
	//	s->len    = 0;
	//	s->tailOp = 0;

	s->size = entries;

	return s;
	}

//----------
//
// edit_script_make_room--
//	Make sure an alignment edit script has enough unused entries available, and
//	enlarge if it doesn't.
//
//----------
//
// Arguments:
//	editscript**	s:			(pointer to) The script to check/enlarge.  If
//								.. reallocation is required we may alter this.
//	u32				entries:	The number of unused entries required.
//	(none)
//
// Returns:
//	Nothing.  Note that a reallocation failure is a fatal error.
//
//----------

static void edit_script_make_room
   (editscript**	_s,
	u32				entries)
	{
	editscript*		s = *_s;

	// do we have enough space already?

	entries += s->len;
	if (s->size >= entries)
		return;  // (yes)

	// reallocate

	entries += entries/2;	// (anticipate 50% future growth)

	*_s = s = realloc_or_die ("edit_script_has_room",
	                          s, edit_script_bytes(entries));

	s->size = entries;
	}

//----------
//
// edit_script_add--
//	Add an repeated operation to an alignment edit script, merging it with the
//	tail of the script if possible.
//
//----------
//
// Arguments:
//	editscript** s:		The script to add to.  If the script has to be
//						.. enlarged, this value may change upon return.
//	u32			 op:	The operation to add.
//	unspos		 rpt:	The repeat count (i.e. how many copies of op to add).
//
// Returns:
//	(nothing).
//
//----------

void edit_script_add
   (editscript**	_s,
	u32				op,
	unspos			rpt)
	{
	editscript*		s = *_s;
	editop*			tail;
	u32				tailRpt;

	// if this operation matches the one currently in the tail, increase the
	// repeat count on the tail;  if the repeat count doesn't have enough count
	// left, fall through to the loop

	if (edit_op_operation(s->tailOp) == op)
		{
		tail    = s->op + s->len-1;
		tailRpt = edit_op_repeat (*tail);

		if (tailRpt + rpt <= maxEditopRepeat)
			{
			*tail = edit_op_add_repeat (*tail, rpt);
			return;
			}
		else
			{
			*tail = edit_op (op, maxEditopRepeat);
			rpt = tailRpt + rpt - maxEditopRepeat;
			}
		}

	// loop, adding new operation(s) to the end of the script

	while (rpt > maxEditopRepeat)
		{
		edit_script_put (_s, op, maxEditopRepeat);
		rpt -= maxEditopRepeat;
		}

	edit_script_put (_s, op, rpt);
	}

//----------
//
// edit_script_put--
//	Add an repeated operation to an alignment edit script, tacking a new entry
//	onto the script.
//
//----------
//
// Arguments:
//	editscript** s:			The script to add to.  If the script has to be
//							.. enlarged, this value may change upon return.
//	u32			 op, rpt:	The operation to add, including a repeat count.
//
// Returns:
//	(nothing).
//
//----------

static void edit_script_put
   (editscript**	_s,
	u32				op,
	u32				rpt)
	{
	editscript*		s;

	edit_script_make_room (_s, 1); // (make sure we have room for one operation)

	s = *_s;
	s->op[s->len++] = edit_op (op, rpt);
	s->tailOp = op;
	}

//----------
//
// edit_script_append--
//	Copy one alignment edit script to the end of another.
//
//----------
//
// Arguments:
//	editscript**	dst:	(pointer to) The script to copy to.  If
//							.. reallocation is required we may alter this.
//	editscript*		src:	(pointer to) The script to copy from.
//
// Returns:
//	A pointer to an empty edit script.  The caller is responsible for disposing
//	of this memory, for which purpose free() can be used.
//
//----------

void edit_script_append
   (editscript**	_dst,
	editscript*		src)
	{
	editscript*		dst;
	editop*			s, *d;
	u32				toCopy;
	u32				sOp;
	u32				sRpt, dRpt;

	if (src->len == 0) return;

	// make sure we have enough room

	edit_script_make_room (_dst, src->len);
	dst = *_dst;

	// copy dst to src

	s  = src->op;
	d  = dst->op + dst->len-1;
	toCopy = src->len;

	sOp = edit_op_operation (*s);
	if (sOp == dst->tailOp)
		{
		dRpt = edit_op_repeat (*d);
		sRpt = edit_op_repeat (*s);
		if (dRpt + sRpt <= maxEditopRepeat)
			*d = edit_op_add_repeat (*d, sRpt);
		else
			{
			*(d++) = edit_op (sOp, maxEditopRepeat);
			*d     = edit_op (sOp, dRpt + sRpt - maxEditopRepeat);
			dst->len++;
			}
		s++;  toCopy--;
		}
	d++;

	memcpy (d, s, toCopy*sizeof(editop));

	dst->len    += toCopy;
	dst->tailOp =  src->tailOp;
	}

//----------
//
// edit_script_reverse--
//	Reverse the items in an alignment edit script, in place.
//
//----------
//
// Arguments:
//	editscript*	s:	The script to modify.
//
// Returns:
//	(nothing)
//
//----------

void edit_script_reverse
   (editscript*	s)
	{
	int			i, j;
	editop		t;

	for (i=0,j=s->len-1 ; i<j ; i++,j--)
		{ t = s->op[i];  s->op[i] = s->op[j];  s->op[j] = t; }
	}

//----------
//
// edit_script_run_of_subs, edit_script_run_of_subs_match--
//	Find the length of the current run of substitutions in an alignment edit
//	script.
//
//----------
//
// Arguments:
//	editscript*	s:		The script being parsed.
//	u32*		opIx:	Current parse location in the script.  Upon return this
//						.. is updated to point to the next operation beyond the
//						.. run.
//	const u8*	p, q:	Pointer to the sequences' nucleotides (corresponding to
//						.. the parse location).  These are used only to provide
//						.. a match count, and can be NULL if match is NULL.
//	unspos*		match:	Place to return the number of nucleotide matches in the
//						.. run (including upper/lower mismatches).
//
// Returns:
//	The length of the run.
//
//----------

u32 edit_script_run_of_subs
   (editscript*	s,
	u32*		_opIx)
	{
	u32			opIx = (u32) *_opIx;
	u32			rpt, run;

	run = 0;
	while ((opIx < s->len) && (edit_op_operation(s->op[opIx]) == editopSub))
		{
		rpt = edit_op_repeat(s->op[opIx]);  opIx++;
		run += rpt;
		}

	*_opIx = opIx;
	return run;
	}


u32 edit_script_run_of_subs_match
   (editscript*	s,
	u32*		_opIx,
	const u8*	p,
	const u8*	q,
	unspos*		_match)
	{
	u32			opIx  = *_opIx;
	unspos		match = *_match;
	u32			rpt, run;
	u8			pCh, qCh;

	run = 0;
	match = 0;
	while ((opIx < s->len)
	    && (edit_op_operation(s->op[opIx]) == editopSub))
		{
		rpt = edit_op_repeat(s->op[opIx]);  opIx++;
		run += rpt;
		while (rpt-- > 0)
			{
            pCh = *(p++);  qCh = *(q++);
            if (dna_toupper(pCh) == dna_toupper(qCh)) match++;
			}
		}

	*_opIx  = opIx;
	*_match = match;
	return run;
	}

//----------
//
// edit_script_indel_len--
//	Find the length of the current "run" of indels in an alignment edit script.
//
//----------
//
// Arguments:
//	editscript*	s:		The script being parsed.
//	u32*		opIx:	Current parse location in the script.  Upon return this
//						.. is updated to point to the next operation beyond the
//						.. indel.
//	unspos*		i, j:	Current parse location in the sequences.  Upon return
//						.. these are updated.
//
// Returns:
//	The length of the run.
//
//----------

u32 edit_script_indel_len
   (editscript*	s,
	u32*		opIx,
	unspos*		i,
	unspos*		j)
	{
	editop		op;
	u32			rpt;

	if (s->len <= (u32) *opIx)
		return 0;

	op = s->op[*opIx];
	rpt = edit_op_repeat(op);

	switch (edit_op_operation(op))
		{
		case editopIns: *j += rpt;  break;
		case editopDel: *i += rpt;  break;
		}

	(*opIx)++;
	return rpt;
	}

//----------
//
// dump_edit_script--
//	Print the raw contents of an alignment edit script.
//
//----------
//
// Arguments:
//	FILE*		f:	The file to print to.
//	editscript* s:	The script to print.
//
// Returns:
//	(nothing).
//
//----------

void dump_edit_script
   (FILE*		f,
	editscript*	s)
	{
	char*		opName[4] = { "???", "INS", "DEL", "SUB" };
	u32			op, rpt;
	u32			ix;

	for (ix=0 ; ix<s->len ; ix++)
		{
		op  = edit_op_operation (s->op[ix]);
		rpt = edit_op_repeat    (s->op[ix]);
		fprintf (f, "%dx%s\n", rpt, opName[op]);
		}
	}

