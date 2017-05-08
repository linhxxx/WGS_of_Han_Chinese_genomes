//#define _EXAM_ASSERT_TEST_
//#define _VERY_GOSSIP_

#include"general.h"
#include"main.h"
#include"func.h"
#include"func2.h"
#include"pileup.h"
#include<iostream>
#include<boost/lexical_cast.hpp>

using namespace std;
using namespace boost;
using namespace aqua;

int Pileup(int mode, MSMUS & indel, MSMUS & snp)
{
	int count(0);
	char isam_cmd[BUFF_LEN];

	for(uint i(0); i < bam_path.size(); ++i)
	{
		cerr<<"Piling up "<<bam_path[i]<<endl;
		snprintf(isam_cmd, BUFF_LEN, "%s view %s 2>/dev/null", samtools_path, bam_path[i].c_str());

		DBG(isam_cmd);
		FILE * ISAM = popen(isam_cmd, "r");

		char buff[BUFF_LEN];
		string sam_content;
		VSTR line_vec;
		vector<pair<uint, char> > cigar_vec;

		while(true)
		{
			fgets(buff, (BUFF_LEN - 1), ISAM);
			if(feof(ISAM))
				break;
			if(strlen(buff) == 0)
				continue;
			sam_content = buff;
			line_vec.clear();
			if(!regex_split(back_inserter(line_vec), sam_content))
			{
				DBG(string("Split Error: ") + sam_content);
				continue;
			}

			ASSERT((!line_vec[0].empty()), "The head column shouldn't be empty!");

			if(line_vec[0][0] == '@')
				continue;
			else if(!(atoi(line_vec[flag].c_str()) & SAM_UNMAPPED))
			{
				//u16_t rd_length(line_vec[seq].size());
				u32_t rd_pos(lexical_cast<u32_t>(line_vec[pos]));

				if(mode & PILEUP_MODE_INDEL)
				{
					cigar_vec.clear();
					__SeperateCigarPure(line_vec[cigar], cigar_vec);

					if(cigar_vec.size() != 1)
					{
						uint accumulation(0);
						uint limit(cigar_vec.size() - 1);
						for(uint i(0); i < limit; ++i)
						{
							u32_t real_left_coordination_inc(rd_pos + accumulation);
							//u32_t real_right_coordination_exc(real_left_coordination_inc + cigar_vec[i].first);
							if(cigar_vec[i].second == STATE_D)
							{
								string & str(indel[line_vec[rname]][real_left_coordination_inc]);
								str += STATE_D;
								str += lexical_cast<string>(cigar_vec[i].first);
								str += '/';
								++count;
							}
							else if(cigar_vec[i].second == STATE_I)
							{
								string & str(indel[line_vec[rname]][real_left_coordination_inc]);
								str += STATE_I;
								str += line_vec[seq].substr(accumulation, cigar_vec[i].first);
								str += '/';
								++count;
								accumulation += cigar_vec[i].first;
							}
							else
							{
								accumulation += cigar_vec[i].first;
								continue;
							}
						}
					}
				}

				if(mode & PILEUP_MODE_SNP)
				{
					ASSERT(line_vec.size()>=12, "Column should larger than 12 when aligned in [Pileup]");

					uint clipAddon(0);
					{
						uint limit(line_vec[cigar].size());
						string pos;
						for(register int j(0); j < limit; ++j)
						{
							if(isalpha(line_vec[cigar][j]))
							{
								if(line_vec[cigar][j] == 'S')
									clipAddon += atoi(pos.c_str());
								break;
							}
							else
								pos += line_vec[cigar][j];
						}
					}
					VEC_UINT_CHAR cigarVec;
					{
						uint j(line_vec.size() - 1);
						string cigarStr;
						for(; j > (tag - 1); --j)
						{
							ASSERT(line_vec[j].size()>2, "Tags should have more than 2 characters in [Pileup]");

							if(line_vec[j][0] == 'M' && line_vec[j][1] == 'D')
							{
								cigarStr = line_vec[j].substr(line_vec[j].rfind(':'), string::npos);
								break;
							}
						}
						__SeperateCigarMismatch(cigarStr, cigarVec);
					}

					{
						uint limit(cigarVec.size());
						if(limit)
						{
							for(register uint i(0); i < limit; ++i)
							{
								u32_t snp_coordination(rd_pos + clipAddon + cigarVec[i].first);
								snp[line_vec[rname]][snp_coordination] += cigarVec[i].second;
								++count;
							}
						}
					}
				}
			}
		}
		pclose(ISAM);
	}

	return count;
}

int PileupPolish(int mode, MSMUS & indel, MSMUS & snp)
{
	int count(0);
	if(mode & PILEUP_MODE_INDEL)
	{
		for(MSMUS_IT idit(indel.begin()); idit != indel.end(); ++idit)
		{
			for(MUS_IT posit(idit->second.begin()); posit != idit->second.end(); ++posit)
			{
				PileupCigarMergeAndSelect(posit->second);
				++count;
			}
			for(MUS_IT posit(idit->second.begin()); posit != idit->second.end();)
			{
				if(posit->second.empty())
				{
					posit = idit->second.erase(posit);
				}
				else
					++posit;
			}
		}
	}
	if(mode & PILEUP_MODE_SNP)
	{
		for(MSMUS_IT idit(snp.begin()); idit != snp.end(); ++idit)
		{
			for(MUS_IT posit(idit->second.begin()); posit != idit->second.end(); ++posit)
			{
				PileupCureSNP(posit->second);
				++count;
			}
			for(MUS_IT posit(idit->second.begin()); posit != idit->second.end();)
			{
				if(posit->second.empty())
				{
					posit = idit->second.erase(posit);
				}
				else
					++posit;
			}
		}
	}

	return count;
}

int PileupOutput(int mode, MSMUS & indel, MSMUS & snp, char * fn)
{
	int count(0);

	fstream O;
	OpenFile(O, fn, "wc-");

	MSMUS_IT idit;
	MSMUS_IT idit_end;
	if(mode & PILEUP_MODE_INDEL)
	{
		idit = indel.begin();
		idit_end = indel.end();
	}
	else if(mode & PILEUP_MODE_SNP)
	{
		idit = snp.begin();
		idit_end = snp.end();
	}
	for(; idit != idit_end; ++idit)
	{
		MUPSS mapDetailBase;
		MSMUS_IT indel_it(indel.find(idit->first));
		if(indel_it != indel.end())
		{
			for(MUS_IT posit(indel_it->second.begin()); posit != indel_it->second.end(); ++posit)
			{
				mapDetailBase[posit->first].first = posit->second;
			}
		}
		MSMUS_IT snp_it(snp.find(idit->first));
		if(snp_it != snp.end())
		{
			for(MUS_IT posit(snp_it->second.begin()); posit != snp_it->second.end(); ++posit)
			{
				mapDetailBase[posit->first].second = posit->second;
			}
		}

		if(mapDetailBase.size())
			O<<">"<<idit->first<<endl;
		for(MUPSS_IT it(mapDetailBase.begin()); it != mapDetailBase.end(); ++it)
		{
			O<<it->first<<'\t'<<(it->second.first.empty()?"-":it->second.first)<<'\t'<<(it->second.second.empty()?"-":it->second.second)<<endl;
			++count;
		}
	}

	return count;
}

inline void PileupCigarMergeAndSelect(string & str)
{
	ASSERT(str.size(), "Input cigar string shouldn't be empty in [PileupCigarMergeAndSelect]");
	typedef unordered_map<string, uint> MAP_STR_UINT;
	MAP_STR_UINT cigar_map;
	string cigar;
	cigar += str[0];
	uint limit(str.size());
	for(register uint i(1); i < limit; ++i)
	{
		if(str[i] == '/')
		{
			++cigar_map[cigar];
			cigar.clear();
			continue;
		}
		cigar += str[i];
	}
	++cigar_map[cigar];

	str.clear();
	MAP_STR_UINT::iterator itLimit(cigar_map.end());

	uint count(0);
	for(MAP_STR_UINT::iterator it(cigar_map.begin()); it != itLimit; ++it)
	{
		//cerr<<"[CigarMergeAndSelect]: Cigar: "<<it->first<<" Count: "<<it->second<<endl;
		//if(it->second < 2)
		//	continue;
		if(it->second > count && it->second >= SUPPORTING_READS_LOW)
		{
			str = it->first;
			str += "/";
			str += lexical_cast<string>(it->second);
			count = it->second;
		}
		else if(it->second == count)
		{
			if(str[0] == STATE_I && it->first[0] == STATE_D)
			{
				str = it->first;
				str += "/";
				str += lexical_cast<string>(it->second);
			}
		}
	}
	//cerr<<"[CigarMergeAndSelect]: FinalCigar: "<<str<<endl;
}

inline void PileupCureSNP(string & str)
{
	typedef unordered_map<char, uint> MAP_CHAR_UINT;
	typedef map<uint, char> MAP_UINT_CHAR;

	double totalAlleleCount(0);

	MAP_CHAR_UINT mapSnp;
	uint limit(str.size());
	for(register uint i(0); i < limit; ++i)
	{
		++mapSnp[str[i]];
	}
	MAP_UINT_CHAR mapSnpSort;
	for(MAP_CHAR_UINT::iterator it(mapSnp.begin()); it != mapSnp.end(); ++it)
	{
		mapSnpSort[it->second] = it->first;
		++totalAlleleCount;
	}
	str.clear();
	try {
		for(MAP_UINT_CHAR::reverse_iterator it(mapSnpSort.rbegin()); it != mapSnpSort.rend(); ++it)
		{
			if(it->first >= SNP_LOWER_LIMIT && it->first / totalAlleleCount >= 0.7)
			{
				str += it->second;
				str += lexical_cast<string>(it->first);
			}
		}
	}
	catch (exception_detail::clone_impl<exception_detail::error_info_injector<bad_lexical_cast> > & e) {
		cerr<<"Bad lexical cast in CureSnp: "<<endl
			<<e.what()<<endl
			<<"Str: "<<str<<endl;
	}
}
