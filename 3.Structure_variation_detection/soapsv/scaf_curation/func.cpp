//#define _EXAM_ASSERT_TEST_
//#define _VERY_GOSSIP_

#include"main.h"
#include"func.h"
#include"general.h"
#include<iostream>
#include<string>
#include<vector>
#include<boost/regex.hpp>
#include<boost/unordered_map.hpp>
#include<map>

using namespace std;
using namespace boost;
using namespace aqua;

int ImportSamList(char * filename, VSTR & filelist_vec)
{
	cerr<<"Importing SAM files list... ";
	int rtValue(ImportToVector(filename, filelist_vec, "\n", true));
	cerr<<"End!"<<endl;
	return rtValue;
}

int ImportSam(char * filename, VSTR & sam_vec)
{
	return ImportToVector(filename, sam_vec, "\n", false);
}

ostream & operator << (ostream & O, const SV_VALIDATION_PARAM_STRUCT & sv)
{
	O<<sv.seq<<'\t'<<sv.begin<<'\t'<<sv.end<<'\t'<<sv.flag<<'\t'<<sv.type<<'\t'<<sv.genotype<<"\n";
	return O;
}

int ProcessSamList(VSTR & filelist_vec, int cpu, char * filename)
{
	fstream O;
	OpenFile(O, filename, "wc-");
	if(!O)
	{
		cerr<<"Output channel blocked!"<<endl;
		exit(EXIT_FAILURE);
	}

	cerr<<"Building vector for multiple threads...";
	vector<SAM_PARAM_STRUCT> param_vector(filelist_vec.size());
	CThreadManager threadMgr(cpu);
	threadMgr.SetTimer(0.001);
	for(register uint i(0); i < filelist_vec.size(); ++i)
	{
		param_vector[i].filename = filelist_vec[i];
		threadMgr.AddThread(_ProcessSam, (void *)(&param_vector[i]));
	}

	filelist_vec.clear();
	cerr<<" End!"<<endl;
	cerr<<param_vector.size()<<" tasks added."<<endl
		<<"Processing..."<<endl;

	threadMgr.Run();

	cerr<<"Starting output..."<<endl;

	for(uint i(0); i < param_vector.size(); ++i)
	{
		O<<param_vector[i].sv;
	}

	cerr<<"Check it out, babe!"<<endl;

	return param_vector.size();
}

inline int __GetSVInfomation(const VSTR & line_sn, SV_VALIDATION_PARAM_STRUCT & sv)
{
	if(line_sn.size() >= 2 && line_sn[1].find("SN:") != string::npos)
	{
		//SN:chr16_21733228_21733233_150_D_5_AGGAA
		RegexSearch match(line_sn[1], "SN:(\\S+?)_(\\d+)_(\\d+)_(\\d+)_(\\w)_(\\d+)_(\\w+)");

		if(!match.is_success())
			return 0;

		sv.seq       = match[1].str();
		sv.begin     = lexical_cast<u32_t>(match[2].str());
		sv.end       = lexical_cast<u32_t>(match[3].str());
		sv.start_pos = lexical_cast<u16_t>(match[4].str());
		sv.type      = match[5].str().empty() ? 'N' : match[5].str()[0];
		sv.genotype  = match[7].str();

		return 1;
	}
	else
		return 0;
}

//Returns the coordination in vector
inline uint __SeperateCigar(string & str, vector<pair<uint, char> > & vec)
{
	string pos;
	uint accumulate(0);
	uint limit(str.size());
	for(register uint i(0); i < limit; ++i)
	{
		if(isalpha(str[i]) || str[i] == '-')
		{
			vec.push_back(make_pair((accumulate), str[i]));
			ASSERT(pos.size());
			accumulate += atoi(pos.c_str());
			pos.clear();
		}
		else
			pos += str[i];
	}

	return vec.size();
}

double __CalcComplexity(const string & str)
{
	if(str.size() < 2)
		return 1;
	u8_t possibility(str.size() - 1);
	if(possibility > 16)
		possibility = 16;
	unordered_map<string, u8_t> map_str;
	uint limit(str.size() - 1);
	for(uint i(0); i < limit; ++i)
	{
		try {
			++map_str[str.substr(i, 2)];
		}
		catch(out_of_range & e) {
			cerr<<"Out of range in __CalcComplexity: "<<endl
				<<e.what()<<endl
				<<"str: "<<str<<endl;
		}
	}

	return ((double)map_str.size()/possibility);
}

void * _ProcessSam(void * param)
{
	static uint process_count(1);
	if(!(process_count % 10000))
		cerr<<"."<<endl;

	SAM_PARAM_STRUCT * param_ptr = (SAM_PARAM_STRUCT *)(param);

	fstream I;
	VSTR sam_content;
	char * fn = strdup(param_ptr->filename.c_str());
	if(!ImportSam(fn, sam_content))
	{
		cerr<<param_ptr->filename<<" missing!"<<endl;
		return (void *)NULL;
	}

	u32_t sv_left_inc(0);
	u32_t sv_right_exc(0);

	u16_t aligned_read_count(0);
	u16_t gapped_read_count(0);
	u16_t contradict_read_count(0);
	u16_t support_linkage_count(0);

	vector<pair<u32_t, u32_t> > vec_coordination;
	map<uint, uint> insertionCoverageMap;

	double complexity_total(0);
	int aside_influence(0);
	bool marker_pending(true);

	for(uint i(0); i < sam_content.size(); ++i)
	{
		VSTR line_vec;
		if(!regex_split(back_inserter(line_vec), sam_content[i]))
		{
			DBG(string("Can't split line: ") + sam_content[i]);
			continue;
		}

		GOSSIP("Marker 1");

		if(marker_pending && (line_vec[0].find("@SQ") != string::npos))
		{
			marker_pending = false;
			GOSSIP("SQ specific");
			if(!__GetSVInfomation(line_vec, param_ptr->sv))
			{
				DBG(string("@SQ header missing in file: ") + param_ptr->filename);
				continue;
			}
			sv_left_inc = (param_ptr->sv.start_pos + 1);
			sv_right_exc = (sv_left_inc + param_ptr->sv.genotype.size());
		}
		else if(!(atoi(line_vec[flag].c_str()) & SAM_UNMAPPED))
		{
			u16_t rd_length(line_vec[seq].size());
			u32_t rd_pos(lexical_cast<u32_t>(line_vec[pos]));

			if((rd_pos + rd_length <= sv_left_inc) || (rd_pos >= sv_right_exc))
			{
				continue;
			}
			GOSSIP("Marker 2");
			if(  !	(\
					 ( (rd_pos < sv_right_exc) && ((rd_pos + AMB_EDGE) >= sv_right_exc) )\
					 ||\
					 ( ((rd_pos + rd_length) > sv_left_inc) && ((rd_pos + rd_length - AMB_EDGE) <= sv_left_inc))
					)
			  )
			{
				++support_linkage_count;
				if(param_ptr->sv.genotype.size() >= 5)
					complexity_total += __CalcComplexity(param_ptr->sv.genotype);
				else
					complexity_total += 1;

				vector<pair<uint, char> > cigar_vec;
				__SeperateCigar(line_vec[cigar], cigar_vec);
				GOSSIP("Marker 3");
				if(cigar_vec.size() == 1)
				{
					if(cigar_vec[i].second == 'D')
						++aligned_read_count;
					if(param_ptr->sv.type == 'I')
					{
						uint limit(rd_pos + rd_length);
						for(register uint i(rd_pos); i < limit; ++i)
							++insertionCoverageMap[i];
						//cerr<<"[_ProcessSam]: Insertion mask between "<<rd_pos<<" and "<<rd_pos+rd_length<<endl;
					}
				}
				else
				{
					for(register uint i(0); i < (cigar_vec.size() - 1); ++i)
					{
						if(cigar_vec[i].second == 'D')
						{
							u32_t left_coordination_inc(lexical_cast<u32_t>(line_vec[pos]) + cigar_vec[i].first);
							u32_t right_coordination_exc(lexical_cast<u32_t>(line_vec[pos]) + cigar_vec[i+1].first);
							GOSSIP("Marker 4");

							//cerr<<"[_ProcessSam]: Gap found @ left: "<<left_coordination_inc<<" beside "<<sv_left_inc<<" Right: "<<right_coordination_exc<<" beside "<<sv_right_exc<<endl;

							if(left_coordination_inc == sv_left_inc && right_coordination_exc == sv_right_exc)
								++gapped_read_count;
							else if((left_coordination_inc < sv_left_inc) && (left_coordination_inc >= (sv_left_inc - 10)))
								vec_coordination.push_back(make_pair(left_coordination_inc, right_coordination_exc));
							else if((left_coordination_inc >= sv_right_exc) && (left_coordination_inc < (sv_right_exc + 10)))
								vec_coordination.push_back(make_pair(left_coordination_inc, right_coordination_exc));
							else
								NULL;
						}
					}
				}
			}
		}
	}

	{
		//Insertion Judgment
		int totalScore(0);
		bool preGap(false);
		for(uint i(sv_left_inc); i < sv_right_exc; ++i)
		{
			if(insertionCoverageMap.find(i) != insertionCoverageMap.end())
			{
				preGap = false;
				totalScore += 1;
			}
			else
			{
				if(preGap)
				{
					totalScore -= 3;
				}
				else
				{
					preGap = true;
					totalScore -= 40;
				}
			}
		}
		//cerr<<"[_ProcessSam]: TotalScore of sv @ "<<param_ptr->sv.begin<<": "<<totalScore<<endl;
		if(totalScore > 0)
			aligned_read_count += insertionCoverageMap[sv_left_inc];
	}

	//Biased_coordination
	if(vec_coordination.size())
	{
		for(;;)
		{
			if(vec_coordination.size() == 1)
			{
				++aside_influence;
				break;
			}

			typedef map<u32_t, uint> MU;
			typedef MU::iterator MUIT;
			typedef map<u32_t, MU> MMU;
			typedef MMU::iterator MMUIT;
			MMU map_coordination;
			uint disAgree(0);

			for(register uint i(0); i < vec_coordination.size(); ++i)
			{
				++map_coordination[vec_coordination[i].first][vec_coordination[i].second];
			}
			for(MMUIT mmuit(map_coordination.begin()); mmuit != map_coordination.end(); ++mmuit)
			{
				disAgree += mmuit->second.size();
			}

			gapped_read_count += vec_coordination.size();

			if(disAgree * 100 / vec_coordination.size() > 50 /*30*/)
			{
				contradict_read_count += vec_coordination.size();
				break;
			}

			break;
		}
	}


	//Low Complexity
	if(complexity_total / support_linkage_count < COMPLEXITY_THRESHOLD)
		param_ptr->sv.flag |= SV_LC;
	GOSSIP("Marker 5");
	//Aside influence
	aside_influence ? param_ptr->sv.flag |= SV_AI : NULL;

	//Miss linkage
	if(!support_linkage_count)
		param_ptr->sv.flag |= SV_ML;
	GOSSIP("Marker 6");
	//Support
	for(;;)
	{
		if((aligned_read_count + gapped_read_count) == 0)
		{
			param_ptr->sv.flag |= SV_F;
			break;
		}
		u16_t * sr(0);
		u16_t * nsr(0);
		sr = &aligned_read_count;
		nsr = &gapped_read_count;
		GOSSIP("Marker 7");
		if(param_ptr->sv.type == 'I')
			NULL;
		else if(param_ptr->sv.type == 'D')
			swap(sr, nsr);
		else
			NR("NR point 2");
		if(*sr < SUPPORTING_READS_LOW)
			param_ptr->sv.flag |= SV_F;
		if(*sr >= SUPPORTING_READS_LOW && *sr < SUPPORTING_READS_HIGH &&  ((*sr * 100 / (*sr + *nsr)) > 50 /*70*/))
			param_ptr->sv.flag |= SV_LS;
		else if(*sr >= SUPPORTING_READS_HIGH &&  ((*sr * 100 / (*sr + *nsr)) > 50 /*70*/))
			param_ptr->sv.flag |= SV_HS;
		else
		{
			param_ptr->sv.flag |= SV_F;
			param_ptr->sv.flag |= SV_HET;
		}
		if(contradict_read_count * 100 / (contradict_read_count + aligned_read_count + gapped_read_count) > 50)
		{
			param_ptr->sv.flag &= (~SV_HS & ~SV_LS);
			param_ptr->sv.flag |= SV_F;
			param_ptr->sv.flag |= SV_BIAS;
		}
		GOSSIP("Marker 8");
		break;
	}
	free(fn);

	return (void*) NULL;
}
