#ifndef FUNC2_H_AQUA_
#define FUNC2_H_AQUA_

#include"general.h"
#include"func.h"
#include"main.h"
#include<string>
#include<vector>
#include<map>
#include<boost/unordered_map.hpp>
#include<pthread.h>

using namespace std;
using namespace boost;

#define halfInterval (readLength - (2 * AMB_EDGE))
#define interval (halfInterval * 2)
#define BUFF_LEN 16384
#define VEC_INIT_SIZE 512
#define COL_LEN 100
#define GAP_LIMIT 0x3F
#define DETERMINE_CONFUSION_STEP 10

const char STATE_I = 'I';
const char STATE_D = 'D';

struct FASTA_BP_ADDITION
{
	string info;
};

struct FASTA_SNP_ADDITION
{
	string info;
};

struct FASTA_BP
{
	FASTA_BP() : ch(0), app(0), snp(0)
	{
		if(pthread_mutex_init(&mutex, NULL) != 0)
		{
			cerr<<"Error initializting mutex lock in FASTA_BP, program terminated."<<endl;
			exit(EXIT_FAILURE);
		}
	}
	FASTA_BP(char ch_t) : ch(ch_t), app(0), snp(0)
	{
		if(pthread_mutex_init(&mutex, NULL) != 0)
		{
			cerr<<"Error initializting mutex lock in FASTA_BP, program terminated."<<endl;
			exit(EXIT_FAILURE);
		}
	}
	~FASTA_BP()
	{
		if(app)
			free(app);
		if(snp)
			free(snp);
	}

	void freeApp()
	{
		pthread_mutex_lock(&mutex);
		if(app)
			free(app);
		app = 0;
		pthread_mutex_unlock(&mutex);
	}

	void freeSnp()
	{
		pthread_mutex_lock(&mutex);
		if(snp)
			free(snp);
		snp = 0;
		pthread_mutex_unlock(&mutex);
	}

	FASTA_BP_ADDITION * InitApp()
	{
		if(!app)
			app = new FASTA_BP_ADDITION;
		return app;
	}

	FASTA_SNP_ADDITION * InitSnp()
	{
		if(!snp)
			snp = new FASTA_SNP_ADDITION;
		return snp;
	}

	void MutexOn()
	{
		pthread_mutex_lock(&mutex);
	}

	void MutexOff()
	{
		pthread_mutex_unlock(&mutex);
	}

	char ch;
	FASTA_BP_ADDITION * app;
	FASTA_SNP_ADDITION * snp;
	pthread_mutex_t mutex;
};

typedef vector<FASTA_BP> VEC_FASTA_BP;
typedef unordered_map<string, VEC_FASTA_BP> MAP_VEC_FASTA_BP;
typedef MAP_VEC_FASTA_BP::iterator MAP_VEC_FASTA_BP_ITER;
typedef vector<uint> VEC_UINT;
typedef unordered_map<string, VEC_UINT> MAP_VEC_UINT;
typedef MAP_VEC_UINT::iterator MAP_VEC_UINT_ITER;
typedef map<uint, string> MAP_UINT_STR;
typedef MAP_UINT_STR::iterator MAP_UINT_STR_ITER;
typedef vector<pair<uint, char> > VEC_UINT_CHAR;
typedef VEC_UINT_CHAR::iterator VEC_UINT_CHAR_ITER;

struct SEARCH_ERROR_PARAM
{
	explicit SEARCH_ERROR_PARAM(MAP_VEC_FASTA_BP & d, MAP_UINT_STR & s):\
				fasta_fn(0), fastq_fn(0), sam_fn(0), data(d), mapSam(s), len(0), fasta_bad_bit(0), fastq_bad_bit(0) {}
	~SEARCH_ERROR_PARAM()
	{
		if(fasta_fn)
			free(fasta_fn);
		if(fastq_fn)
			free(fastq_fn);
		if(sam_fn)
			free(sam_fn);
	}

	void set(const string & n, uint l)
	{
		name = n; locus = l;
	}

	string name;
	char * fasta_fn;
	char * fastq_fn;
	char * sam_fn;
	MAP_VEC_FASTA_BP & data;
	MAP_UINT_STR & mapSam;
	uint locus;
	u8_t len:6;
	u8_t fasta_bad_bit:1;
	u8_t fastq_bad_bit:1;


};

int ImportFasta(char *, MAP_VEC_FASTA_BP &, bool);

int ExportFasta(char *, MAP_VEC_FASTA_BP &, bool);

extern int _ExportFastaSeq(MAP_VEC_FASTA_BP &, const string &, uint, uint, string &, bool);

void _WrapString(string &, uint);

int ExportDetails(char *, MAP_VEC_FASTA_BP &, bool);

extern int _ExportEachDetails(MAP_VEC_FASTA_BP &, const string &, uint, uint, string &);

int ImportSuspiciousList(char *, MAP_VEC_FASTA_BP & , bool, int, MAP_VEC_UINT &);

int GlobalDetermineSuspicious(MAP_VEC_FASTA_BP & data, int, MAP_VEC_UINT &);

void * _ProcessSamCore(void *);

int ImportSAMOfWholeScaffold(MAP_UINT_STR &, string);

int SearchError(MAP_VEC_FASTA_BP &, MAP_VEC_UINT &, int);

void * _SearchError(void *);

void * __FastaPipe(void *);

void * __FastqPipe(void *);

void * __SamPipe(void *);

void * __DetermineErrorsFromSam(void *);

void CigarMergeAndSelect(string &);

void CureSNP(string &);

int Base2BaseCuration(MAP_VEC_FASTA_BP &, int);

int EliminateError(MAP_VEC_FASTA_BP &);

extern uint __SeperateCigarPure(string &, VEC_UINT_CHAR &);

extern uint __SeperateCigarMismatch(string & str, VEC_UINT_CHAR &);

int SimulateIndel(MAP_VEC_FASTA_BP &);

#endif /* FUNC2_H_AQUA_ */
