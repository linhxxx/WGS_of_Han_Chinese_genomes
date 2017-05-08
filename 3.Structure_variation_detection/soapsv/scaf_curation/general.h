/*
 *      Filename: general.cpp
 *
 *
 *      Description:
 *         Basic functions
 *
 *      Including Marcos:
 *
 *
 *  	Created on: Feb 8, 2010
 *      Author: Ruibang Luo, BGI
 *
 *     	History:
 *         1.
 */

#pragma once
#ifndef GENERAL_H_AQUA_
#define GENERAL_H_AQUA_

#ifndef __cplusplus
#error Please use C++ complier!
#endif

#include<fstream>
#include"gzstream.h"
#include<string>
#include<vector>
#include<boost/lexical_cast.hpp>
#include<boost/regex.hpp>
#include<boost/utility.hpp>
#include<boost/static_assert.hpp>

using namespace std;
using namespace boost;

//Types************************************************************************
typedef vector<string> VSTR;

typedef unsigned int       uint;
typedef unsigned char      uchar;
typedef unsigned short     ushort;
typedef unsigned long      ulong;
typedef unsigned long long ullong;

template <int i> void size_of_integer_equal_to_4() {
  BOOST_STATIC_ASSERT(sizeof(i) == 4);
}

typedef unsigned char      BYTE;
typedef unsigned short     WORD;
typedef unsigned int       DWORD;

typedef unsigned char      u8_t;
typedef unsigned short     u16_t;
typedef unsigned int       u32_t;
typedef unsigned long      u64_t;

//*****************************************************************************

//Debugging********************************************************************
/* Please define _EXAM_ASSERT_TEST_ to enable ASSERTMENT */
#ifdef _VERY_GOSSIP_

	#define GOSSIP(message)\
			{\
				cerr<<(message)<<endl;\
			}

#else

	#define GOSSIP(message) NULL;

#endif /* _VERY_GOSSIP_ */

#ifdef _EXAM_ASSERT_TEST_

	#define DEBUG_INFO (string("[") + __FILE__ + "]" + "[" + __FUNCTION__ + "]" + "[" + lexical_cast<string>(__LINE__) + "]")

	void exam_assert(const char *, const char *, uint);
	void exam_assert(const char *, const char *, uint, const char *);
	#define ASSERT(condition, ...) \
			{\
				if(!(condition))\
					exam_assert(__FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__);\
			}

	#define DBG(message)\
			{\
				cerr<<(message)<<endl;\
			}

	#define NR(message)\
			{\
				cerr<<DEBUG_INFO<<" Will never reach here. "<<(message)<<endl;\
			}

#else

	#define ASSERT(condition, ...) NULL;
	#define DEBUG_INFO NULL;
	#define DBG(message) NULL;
	#define NR NULL;

#endif /* _EXAM_ASSERT_TEST_ */
//IO operators and manipulators************************************************
namespace aqua
{
ostream & line(ostream &);
ostream & endt(ostream &);
}
//*****************************************************************************

//*****************************************************************************

	//File_operation***************************************************************
namespace aqua
{
int _JudgeStdio(char *&, const string &);
_Ios_Openmode _DefineOpenMode(const string &);
int FileExist(const char *, const bool);
/* File open operators*/
/*
+: append
b: binary
r: input
w: output
-: truncate
c: check for existence, exit if error.
*/
int OpenFile(fstream &, char *, const string);
int OpenFileGZ(gzstreambase &, char *, const string);
}
//*****************************************************************************

//Regex************************************************************************
namespace aqua
{
class RegexBasic : boost::noncopyable
{
public:
	explicit RegexBasic(const string & reg_str) : reg(reg_str), success(false) {}

protected:
	boost::regex reg;
	bool success;
	boost::smatch matched;
};

class RegexSearch : public RegexBasic
{
public:
	explicit RegexSearch(const string & target_str, const string reg_str) : RegexBasic(reg_str)
	{
		success = boost::regex_search(target_str, matched, reg);
	}

	bool is_success();

	smatch::const_reference operator [](uint);

	smatch::size_type size();
};
}

namespace aqua
{
int StrToVector(VSTR &, string &, const string &);
int ImportToVector(char *, VSTR &, const string &, const bool);
int ImportToVectorGZ(char *, VSTR &, const string &, const bool);
int SplitByRegex(VSTR &, string &, const char *);
int SplitBySS(VSTR &, string &, const char);
}
//*****************************************************************************

namespace aqua
{
string SuitFilenameToShell(string);
}

#endif /* GENERAL_H_ */
