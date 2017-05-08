#ifndef PILEUP_H_AQUA_
#define PILEUP_H_AQUA_

#include<string>
#include<boost/unordered_map.hpp>
#include<map>
#include<vector>

using namespace std;
using namespace boost;

#define PILEUP_MODE_INDEL	0x1
#define PILEUP_MODE_SNP		0x2

typedef unordered_map<uint, string> UMAP_UINT_STR;
typedef UMAP_UINT_STR MUS;
typedef MUS::iterator MUS_IT;
typedef unordered_map<string, MUS> UMAP_STR_MAP_UINT_STR;
typedef UMAP_STR_MAP_UINT_STR MSMUS;
typedef MSMUS::iterator MSMUS_IT;

typedef map<uint, pair<string, string> > MUPSS;
typedef MUPSS::iterator MUPSS_IT;

int Pileup(int, MSMUS &, MSMUS &);
int PileupPolish(int, MSMUS &, MSMUS &);
int PileupOutput(int, MSMUS &, MSMUS &, char *);
inline void PileupCigarMergeAndSelect(string &);
inline void PileupCureSNP(string &);

#endif /* PILEUP_H_ */
