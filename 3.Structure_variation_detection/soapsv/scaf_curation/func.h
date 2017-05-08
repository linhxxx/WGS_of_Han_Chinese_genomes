#pragma once
#ifndef FUNC_H_AQUA_
#define FUNC_H_AQUA_

#define _EXAM_ASSERT_TEST_
#include<string>
#include<vector>
#include"general.h"
#include"threadmanager.h"

using namespace std;

#define SAM_SELF_UNMAPPED 0x04
#define SAM_MATE_UNMAPPED 0X08
#define SAM_UNMAPPED ((SAM_SELF_UNMAPPED) | (SAM_MATE_UNMAPPED))
#define AMB_EDGE 7
#define COMPLEXITY_THRESHOLD 0.5

enum sam_column {
	qname = 0, flag, rname, pos, mapq, cigar, mrnm, mpos, isize, seq, qual, tag
};

#define SV_F  0x01
#define SV_LS 0x02
#define SV_HS 0x04
#define SV_ML 0x08
#define SV_LC 0x10
#define SV_AI 0x20
#define SV_BIAS 0x40
#define SV_HET 0x80

int ImportSamList(char *, VSTR &);

int ImportSam(char *, VSTR &);

struct SV_VALIDATION_PARAM_STRUCT
{
	SV_VALIDATION_PARAM_STRUCT() : begin(0), end(0), flag(0), type(0), start_pos(0) {}

	string seq;
	u32_t begin;
	u32_t end;
	u16_t flag;
	char type;
	string genotype;
	u16_t start_pos;
};

struct SAM_PARAM_STRUCT
{
	string filename;
	SV_VALIDATION_PARAM_STRUCT sv;
};

ostream & operator << (ostream &, const SAM_PARAM_STRUCT &);

extern int __GetSVInfomation(const VSTR &, SV_VALIDATION_PARAM_STRUCT &);

extern uint __SeperateCigar(string &, vector<pair<uint, char> > &);

extern double __CalcComplexity(const string &);

int ProcessSamList(VSTR &, int, char *);

void * _ProcessSam(void *);

inline void __ProcessSamCore(VSTR &);

#endif /* FUNC_H_AQUA_ */
