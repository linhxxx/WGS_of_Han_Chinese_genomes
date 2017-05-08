//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: edit_script.h
//
//----------

#ifndef edit_script_H			// (prevent multiple inclusion)
#define edit_script_H

// other files

#include "sequences.h"			// sequence stuff

// establish ownership of global variables

#ifdef edit_script_owner
#define global
#else
#define global extern
#endif

//----------
//
// data structures and types
//
//----------

// linked list of alignments

typedef struct alignel
	{
	struct alignel*	next;
	int				isTrivial;
	unspos			beg1, beg2, // (origin-1)
					end1, end2;
	score			s;
	u8*				seq1, *seq2;
	struct editscript* script;
	} alignel;

// edit scripts;  a list of insert, delete and substitute operations

typedef u32 editop;

enum
	{
	editopIns  = 0x1,
	editopDel  = 0x2,
	editopSub  = 0x3
	};

typedef struct editscript
	{
	u32		size;				// the number of entries allocated for op[]
	u32		len;				// the number of entries used
	editop	tailOp;				// most recent operation added
	editop	op[1];				// variable-length array of edit operations
	} editscript;

#define edit_script_ins(s,rpt)  edit_script_add(s,editopIns,rpt)
#define edit_script_del(s,rpt)  edit_script_add(s,editopDel,rpt)
#define edit_script_sub(s,rpt)  edit_script_add(s,editopSub,rpt)

#define edit_script_bytes(entries)	(sizeof(editscript)+(((entries)-1)*sizeof(editop)))

//----------
//
// prototypes for routines in edit_script.c
//
//----------

void        free_align_list               (alignel* a);
editscript* edit_script_new               (void);
void        edit_script_add               (editscript** s, u32 op, unspos rpt);
void        edit_script_append            (editscript** dst, editscript* src);
void        edit_script_reverse           (editscript* s);
u32         edit_script_run_of_subs       (editscript* s, u32* opIx);
u32         edit_script_run_of_subs_match (editscript* s, u32* opIx,
                                           const u8* p, const u8* q, unspos* match);
u32         edit_script_indel_len         (editscript* s,
                                           u32* opIx, unspos* i, unspos* j);
void        dump_edit_script              (FILE* f, editscript* s);

#undef global
#endif // edit_script_H
