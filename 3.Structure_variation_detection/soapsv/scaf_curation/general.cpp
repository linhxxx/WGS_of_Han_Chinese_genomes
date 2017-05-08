/*
 *      Filename: general.cpp
 *
 *
 *      Description:
 *         Basic functions
 *
 *      Including Functions:
 *
 *
 *  	Created on: Feb 8, 2010
 *      Author: Ruibang Luo, BGI
 *
 *     	History:
 *         1.
 */


#define _EXAM_ASSERT_TEST_
#include<iostream>
#include<fstream>
#include"gzstream.h"
#include<string>
#include"general.h"
#include<boost/regex.hpp>
#include<cstring>
#include<sstream>

#define VEC_INIT_SIZE 0xFFU

using namespace std;
using namespace boost;

inline int aqua::FileExist(const char * filename, const bool aggresive)
{
	ifstream I(filename, ios::in | ios::binary);
	if(!I)
	{
		if(aggresive)
		{
			cerr<<endl
				<<DEBUG_INFO
				<<"File missing: "<<filename<<endl;
			exit(EXIT_FAILURE);
		}
		return 0;
	}
	I.close();
	return 1;
}

/* File open operators*/
/*
+: append
b: binary
r: input
w: output
-: truncate
c: check for existence, exit if error.
*/

inline int aqua::_JudgeStdio(char *& file_name, const string & str)
{
	if(!string(file_name).empty() && string(file_name)[0] == '-')
	{
		if(str.find("r") != string::npos)
		{
			cerr<<"Input redirect to STDIN"<<endl;
			file_name = strdup("/dev/stdin");
		}
		else if(str.find("w") != string::npos)
		{
			cerr<<"Output redirect to STDOUT"<<endl;
			file_name = strdup("/dev/stdout");
		}
		else
			return 0;
	}

	return 1;
}

ostream & aqua::line(ostream & O)
{
  return O<<"\n--------------------------------------------------------------------------------\n";
}

ostream & aqua::endt(ostream & O)
{
  return O<<"\t";
}

inline _Ios_Openmode aqua::_DefineOpenMode(const string & str)
{
	_Ios_Openmode open_mode(ios::trunc ^ ios::trunc);

	for(uint i(0); i < str.size(); ++i)
	{
		switch(str[i])
		{
		case '+': open_mode |= ios::app; break;
		case 'b': open_mode |= ios::binary; break;
		case 'r': open_mode |= ios::in; break;
		case 'w': open_mode |= ios::out; break;
		case '-': open_mode |= ios::trunc; break;
		//default: break;
		}
	}
	ASSERT( (!((open_mode & ios::app) && (open_mode & ios::trunc))), "Append and Truncate at the same time.");
	ASSERT( (!((open_mode & ios::in) && (open_mode & ios::out))), "Input and Output at the same time.");
	ASSERT( (!((open_mode & ios::in) && (open_mode & ios::trunc))), "Input and Truncate at the same time.");
	ASSERT( (!((open_mode & ios::in) && (open_mode & ios::app))), "Input and Append at the same time.");

	return open_mode;
}

int aqua::OpenFile(fstream & FH, char * file_name, const string str)
{
	_JudgeStdio(file_name, str);
	FH.open(file_name, _DefineOpenMode(str));
	return FileExist(file_name, (str.find('c') != string::npos));
}

int aqua::OpenFileGZ(gzstreambase & FH, char * file_name, const string str)
{
	_JudgeStdio(file_name, str);
	FH.open(file_name, _DefineOpenMode(str));
	return FileExist(file_name, (str.find('c') != string::npos));
}

int aqua::StrToVector(VSTR & vec_str, string & whole_str, const string & reg_str)
{
	ASSERT((reg_str.size() != 0), "Split Regex equal to 0!");
	char reg = reg_str[0];
	SplitBySS(vec_str, whole_str, reg);

	/*if(vec_str.size() && vec_str[vec_str.size() - 1].empty())
		vec_str.resize(vec_str.size() - 1);*/

	return vec_str.size();
}

int aqua::ImportToVector(char * filename, VSTR & vec_str, const string & reg_str, const bool ambiguous)
{
	fstream I;
	if(!OpenFile(I, filename, (ambiguous ? "rc" : "r")))
		return 0;

	register uint pos(0);
	vec_str.resize(VEC_INIT_SIZE);
	char reg(reg_str[0]);
	while(true)
	{
		getline(I, vec_str[pos], reg);
		if(!I)
		{
			vec_str.resize(pos);
			break;
		}
		++pos;
		if(pos == vec_str.size())
			vec_str.resize(vec_str.size() << 1);
	}

	return vec_str.size();
}

int aqua::ImportToVectorGZ(char * filename, VSTR & vec_str, const string & reg_str, const bool ambiguous)
{
	igzstream I;
	if(!OpenFileGZ(I, filename, (ambiguous ? "rc" : "r")))
		return 0;

	register uint pos(0);
	vec_str.resize(VEC_INIT_SIZE);
	char reg(reg_str[0]);
	while(true)
	{
		getline(I, vec_str[pos], reg);
		if(!I)
		{
			vec_str.resize(pos);
			break;
		}
		++pos;
		if(pos == vec_str.size())
			vec_str.resize(vec_str.size() << 1);
	}

	return vec_str.size();
}

int aqua::SplitByRegex(VSTR & vec, string & str, const char * reg)
{
	regex e(reg);
	regex_split(back_inserter(vec), str, e);
	return vec.size();
}

int aqua::SplitBySS(VSTR & vec, string & str, const char reg)
{

	stringstream ss(str);
	vec.resize(VEC_INIT_SIZE);
	register uint pos(0);
	for(;;)
	{
		getline(ss, vec[pos], reg);
		if(!ss)
		{
			vec.resize(pos);
			break;
		}
		++pos;
		if(pos == vec.size())
			vec.resize(vec.size() << 1);
	}

	return vec.size();
}

smatch::const_reference aqua::RegexSearch::operator [](uint i)
{
	return matched[i];
}

smatch::size_type aqua::RegexSearch::size()
{
	return matched.size();
}

bool aqua::RegexSearch::is_success()
{
	return success;
}

void exam_assert(const char * file_name, const char * func_name, uint line_no)
{
	cerr<<endl<<DEBUG_INFO<<":Assert failure"<<endl;
	abort();
}

void exam_assert(const char * file_name, const char * func_name, uint line_no, const char * reason)
{
	cerr<<endl<<DEBUG_INFO<<":Assert failure\nReason: "<<endl<<reason<<endl;
	abort();
}

string aqua::SuitFilenameToShell(string name)
{
	for(register int i(0); i < name.size();)
	{
		i = name.find_first_of("<>*&$\\;|\'\"[]=:", i);
		if(i != string::npos)
		{
			name.insert(i, "\\");
			i += 2;
		}
	}

	return name;
}
