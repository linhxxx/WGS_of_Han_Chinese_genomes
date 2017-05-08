//#define _EXAM_ASSERT_TEST_
//#define _VERY_GOSSIP_
#include"general.h"
#include"main.h"
#include"func.h"
#include"func2.h"
#include"threadmanager.h"
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<iostream>
#include<sstream>
#include<map>
#include<functional>
#include<boost/regex.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/unordered_map.hpp>
#include<sys/vfs.h>

using namespace std;
using namespace aqua;
using namespace boost;

//SamProcessingSlot
void * (*__SamProcessingSlot)(void *)(NULL);

int ImportFasta(char * file_name, MAP_VEC_FASTA_BP & data, bool ambiguous)
{
	fstream I;
	if(!OpenFile(I, file_name, (ambiguous ? "rc" : "r")))
		return 0;

	string id;
	string seq;
	int count(0);

	getline(I, id, '>');
	while(true)
	{
		getline(I, id, '\n');
		if(!I)
			break;

		getline(I, seq, '>');
		if(seq.empty())
			continue;
		{
			VSTR id_vec;
			int flag1(SplitByRegex(id_vec, id, "\\s+"));
			ASSERT(flag1 != 0, "ID split error in [ImportFasta]");
			id = id_vec[0];
		}
		cerr<<"Importing sequence "<<id<<endl;
		data[id].resize(seq.size());
		++count;
		MAP_VEC_FASTA_BP::iterator it(data.find(id));
		ASSERT((it != data.end()), "Iterator construction failure in [ImportFasta]");

		{
			uint i(0), j(0);
			for(; i < seq.size(); ++i)
			{
				if(seq[i] == '\n' || seq[i] == ' ')
					continue;
				it->second[j].ch = seq[i];
				//cerr<<"[ImportFasta]: "<<j<<'\t'<<it->second[j].ch<<endl;
				++j;
			}
			ASSERT((j <= i));
			it->second.resize(j);
		}

		ASSERT(it->second.size());
		it->second[0].InitApp()->info += (string("S") + lexical_cast<string>(it->second.size()) + " ");
	}

	cerr<<count<<" sequences imported."<<endl;

	return count;
}

int ExportFasta(char * file_name, MAP_VEC_FASTA_BP & data, bool ambiguous)
{
	fstream O;
	OpenFile(O, file_name, (ambiguous ? "wc-" : "w-"));
	int count(0);

	for(MAP_VEC_FASTA_BP_ITER it(data.begin()); it != data.end(); ++it)
	{
		cerr<<"Exporting sequence "<<it->first<<endl;
		++count;
		string seq;
		int flag1;
		flag1 = _ExportFastaSeq(data, it->first, 0, it->second.size(), seq, false);
		ASSERT(flag1, "_ExportFastaSeq returned a empty string!");
		_WrapString(seq, COL_LEN);
		ASSERT(it->second[0].app, "it->second[0].app in [ExportFasta]");
		O<<">"<<it->first<<'\t'<<it->second[0].app->info<<"\n"<<seq;
	}

	return count;
}

inline int _ExportFastaSeq(MAP_VEC_FASTA_BP & data, const string & name, uint start_inc, uint end_exc, string & output_str, bool raw)
{
	output_str.clear();
	MAP_VEC_FASTA_BP_ITER it(data.find(name));
	ASSERT((it != data.end()), "it != data.end() in [_ExportFastaSeq]")
	if(end_exc > it->second.size())
		end_exc = it->second.size();
	DBG(string("Start: ") + lexical_cast<string>(start_inc) + " End: " + lexical_cast<string>(end_exc) + " in [_ExportFastaSeq]");
	ASSERT((start_inc < end_exc), "start_inc < end_exc in [_ExportFastaSeq]");
	for(register uint i(start_inc); i < end_exc; ++i)
	{
		try {
			if(!raw && it->second[i].app && i)
			{
				string * appStr(&(it->second[i].app->info));
				ASSERT(appStr->size(), "FASTA_ADDITION_STRUCT, empty string in [ExportFasta]");
				if((*appStr)[0] == STATE_I)
				{
					string appendStr(appStr->substr(1, appStr->find("/", 1) - 1));
					bool bypass(false);
					for(register uint i(0); i < appendStr.size(); ++i)
					{
						if(isdigit(appendStr[i]))
						{
							cerr<<"Digit in inserted string at: "<<it->first<<'\t'<<i<<'\t'<<*appStr<<endl;
							bypass = true;
							break;
						}
					}
					if(!bypass)
						output_str.append(appendStr);
					//cerr<<(appStr->substr(1, appStr->find("/", 1) - 1) + " inserted at <1-coordinate>" + lexical_cast<string>(i+1))<<endl;
				}
				else if((*appStr)[0] == STATE_D)
				{
					i += (lexical_cast<uint>(appStr->substr(1, appStr->find("/", 1) - 1)));
					//cerr<<(lexical_cast<string>(appStr->substr(1, appStr->find("/", 1) - 1))) + " base bypassed at (1-coordinate)" + lexical_cast<string>(i+1))<<endl;
				}
				output_str += it->second[i].ch;
			}
			else
				output_str += it->second[i].ch;
		}

		catch (exception_detail::clone_impl<exception_detail::error_info_injector<bad_lexical_cast> > & e) {
			cerr<<"Bad lexical cast in _ExportFastaSeq: "<<endl
				<<e.what()<<endl
				<<"it->second[i].app->info: "<<it->second[i].app->info<<endl;
		}
		catch (out_of_range & e) {
			cerr<<"Out of range in _ExportFastaSeq: "<<endl
				<<e.what()<<endl
				<<"it->second[i].app->info: "<<it->second[i].app->info<<endl;
		}
	}

	return output_str.size();
}

void _WrapString(string & str, uint col_len)
{
	string str2;
	uint i(0);
	uint limit(str.size());
	for(; i < limit; i += col_len)
	{
		str2 += str.substr(i, col_len);
		str2 += "\n";
	}
	str = str2;
}

int ExportDetails(char * file_name, MAP_VEC_FASTA_BP & data, bool ambigious)
{
	fstream O;
	OpenFile(O, file_name, (ambigious ? "wc-" : "w-"));
	int count(0);

	for(MAP_VEC_FASTA_BP_ITER it(data.begin()); it != data.end(); ++it)
	{
		cerr<<"Exporting details of "<<it->first<<endl;
		++count;
		string seq;
		int flag1;
		flag1 = _ExportEachDetails(data, it->first, 1, it->second.size(), seq);   // 0 is used for storing the sequence details
		//ASSERT(flag1, "__ExportEachDetails returned a empty string!");
		O<<">"<<it->first<<"\n"<<seq;
	}

	cerr<<"Details of "<<count<<" sequences exported"<<endl;
	return count;
}

inline int _ExportEachDetails(MAP_VEC_FASTA_BP & data, const string & name, uint start_inc, uint end_exc, string & output_str)
{
	output_str.clear();
	MAP_VEC_FASTA_BP_ITER it(data.find(name));
	ASSERT((it != data.end()), "it != data.end() in [_ExportFastaSeq]")
	if(end_exc > it->second.size())
		end_exc = it->second.size();
//cerr<<it->first<<'\t'<<start_inc<<'\t'<<end_exc<<endl;
	ASSERT((start_inc < end_exc));
	for(register uint i(start_inc); i < end_exc; ++i)
	{
		if(it->second[i].app || it->second[i].snp)
		{
			output_str += (it->first + "\t" + lexical_cast<string>(i+1) + "\t");
			if(it->second[i].app)
			{
				string & appStr(it->second[i].app->info);
				ASSERT(appStr.size(), "FASTA_ADDITION_STRUCT, empty string in [ExportFasta]");
				output_str += appStr;
			}
			if(it->second[i].snp)
			{
				string & snpStr(it->second[i].snp->info);
				ASSERT(snpStr.size(), "FASTA_SNP, empty string in [ExportFasta]");
				output_str += "/";
				output_str += snpStr;
			}
			output_str += "\n";
		}
	}

	return output_str.size();
}

int ImportSuspiciousList(char * file_name, MAP_VEC_FASTA_BP & data, bool ambigious, int readLength, MAP_VEC_UINT & suspiciousLoci)
{
	fstream I;
	if(!OpenFile(I, file_name, (ambigious ? "rc" : "r")))
		return 0;

	int errorAccumulate(0);
	int count(0);

	string strtmp;
	while(true)
	{
		if(errorAccumulate > 100)
		{
			cerr<<"Too many errors in the suspicious list file, please check and run again!"<<endl;
			exit(EXIT_FAILURE);
		}

		getline(I, strtmp, '\n');
		if(!I)
			break;
		VSTR vectmp;
		MAP_VEC_FASTA_BP_ITER it;
		SplitByRegex(vectmp, strtmp, "\\s+");
		it = data.find(vectmp[0]);
		uint left(lexical_cast<uint>(vectmp[1]) - 1);
		if((it == data.end()) || (vectmp.size() < 2) || left == 0 || left >= it->second.size())
		{
			++errorAccumulate;
			cerr<<"Suspicious locus: "<<strtmp<<" missing in fasta."<<endl;
			continue;
		}

		if(left < uint(interval))
			left = interval;
		if(vectmp.size() == 2)
			suspiciousLoci[vectmp[0]].push_back(left);
		else if(vectmp.size() == 3)
		{
			uint right(lexical_cast<uint>(vectmp[2]) - 1);
			if(right >= (it->second.size() - interval))
			{
				right = (it->second.size() - interval);
				suspiciousLoci[vectmp[0]].push_back(right);
			}
			if(right < left || right > it->second.size())
			{
				++errorAccumulate;
				cerr<<"Suspicious locus: "<<strtmp<<" missing in fasta."<<endl;
				continue;
			}

			for(; left < right; left += interval)
				suspiciousLoci[vectmp[0]].push_back(left);
		}
		else
		{
			++errorAccumulate;
			continue;
		}

		++count;
	}
	return count;
}

int GlobalDetermineSuspicious(MAP_VEC_FASTA_BP & data, int readLength, MAP_VEC_UINT & suspiciousLoci)
{
	int count(0);
	for(MAP_VEC_FASTA_BP_ITER it(data.begin()); it != data.end(); ++it)
	{
		string details;
		suspiciousLoci[it->first].clear();
		MAP_VEC_UINT_ITER sLit(suspiciousLoci.find(it->first));
		if(it->second.size() < uint(2 * interval))
		{
			details += "Bypassed ";
		}
		else
		{
			uint limit(it->second.size() - interval);
			uint i(interval);
			for( ;i < limit; i += interval)
			{
				if(it->second[i].ch != 'N')
					sLit->second.push_back(i);
				++count;
			}
			if(it->second[it->second.size() - interval].ch != 'N')
				sLit->second.push_back(it->second.size() - interval);
		}

		ASSERT(it->second[0].app, "Fasta [0] app not initialized!");
		it->second[0].app->info += details;
	}

	cerr<<count<<" suspicious loci have been marked."<<endl;

	return count;
}

void * _ProcessSamCore(void * param)
{
	SEARCH_ERROR_PARAM * p = (SEARCH_ERROR_PARAM*) param;

	u32_t sv_left_inc(interval);
	u32_t sv_right_exc(interval + p->len);
	uint zeroCoordination((p->locus - interval));

	//cerr<<"[_ProcessSamCore]: sv_left_inc/sv_right_exc/locus/zeroCoordination: "<<sv_left_inc<<'/'<<sv_right_exc<<'/'<<p->locus<<'/'<<zeroCoordination<<endl;

	u16_t aligned_read_count(0);
	u16_t gapped_read_count(0);
	u16_t contradict_read_count(0);
	u16_t support_linkage_count(0);

	vector<pair<u32_t, u32_t> > vec_coordination;

	double complexity_total(0);
	int aside_influence(0);
	int bias_influence(0);
	string sv_genotype;

	string & info(p->data[p->name][p->locus].app->info);
	ASSERT(info.size(), "App should be already allocated in [_ProcessSamCore]");
	char sv_type(info[0]);
	unordered_map<string, uint> insertionCountingMap;

	u8_t sv_flag(0);

	string tmpStr;
	ifstream I(p->sam_fn);
	while(true)
	{
		getline(I, tmpStr, '\n');
		if(!I)
			break;
		VSTR line_vec;
		if(!regex_split(back_inserter(line_vec), tmpStr))
		{
			DBG(string("Split Error: ") + tmpStr);
			continue;
		}

		//GOSSIP("Marker 1 in [_ProcessSamCore]");
		ASSERT((!line_vec[0].empty()), "The head column shouldn't be empty!");

		if(line_vec[0][0] == '@')
			continue;
		else if(!(atoi(line_vec[flag].c_str()) & SAM_UNMAPPED))
		{
			u16_t rd_length(line_vec[seq].size());
			u32_t rd_pos(lexical_cast<u32_t>(line_vec[pos]) - 1);

			//cerr<<"[_ProcessSamCore]: rd_pos: "<<rd_pos<<endl;

			if((rd_pos + rd_length <= sv_left_inc) || (rd_pos >= sv_right_exc))
			{
				continue;
			}
			//GOSSIP("Marker 2 in [_ProcessSamCore]");
			if(  !	(\
					 ( (rd_pos < sv_right_exc) && ((rd_pos + AMB_EDGE) >= sv_right_exc) )\
					 ||\
					 ( ((rd_pos + rd_length) > sv_left_inc) && ((rd_pos + rd_length - AMB_EDGE) <= sv_left_inc))
					)
			  )
			{
				++support_linkage_count;

				if(sv_type == STATE_I && rd_pos <= sv_left_inc && (rd_pos + rd_length)>= sv_right_exc)
					++insertionCountingMap[line_vec[seq].substr(sv_left_inc - rd_pos, p->len)];
				_ExportFastaSeq(p->data, p->name, zeroCoordination, zeroCoordination + p->len, sv_genotype, true);

				if(p->len >= 5)
				{
					try {
						if(sv_type == STATE_D)
						{
							NULL;
						}
						else if(sv_type == STATE_I && rd_pos <= sv_left_inc && (rd_pos + rd_length)>= sv_right_exc)
						{
							sv_genotype = line_vec[seq].substr(sv_left_inc - rd_pos, p->len);
						}
						else
							complexity_total += 1;
					}
					catch(out_of_range & e) {
						cerr<<"Out of range in _ProcessSamCore: "<<endl
							<<e.what()<<endl
							<<"line_vec[seq]: "<<line_vec[seq]<<endl
							<<"Rd_pos: "<<rd_pos<<endl
							<<"Left_sv_inc: "<<sv_left_inc<<endl
							<<"Right_sv_exc: "<<sv_right_exc<<endl;
					}
					complexity_total += __CalcComplexity(sv_genotype);
				}
				else
					complexity_total += 1;

				vector<pair<uint, char> > cigar_vec;
				__SeperateCigar(line_vec[cigar], cigar_vec);
				//GOSSIP("Marker 3 in [_ProcessSamCore]");
				if(cigar_vec.size() == 1)
					++aligned_read_count;
				else
				{
					uint limit(cigar_vec.size() - 1);
					for(register uint i(0); i < limit; ++i)
					{
						if(cigar_vec[i].second == STATE_D || cigar_vec[i].second == STATE_I)
						{
							u32_t left_coordination_inc(rd_pos + cigar_vec[i].first);
							u32_t right_coordination_exc(rd_pos + cigar_vec[i+1].first);
							//GOSSIP("Marker 4 in [_ProcessSamCore]");

							//cerr<<"[_ProcessSamCore]: Gap found @ left: "<<left_coordination_inc<<" beside "<<sv_left_inc<<" Right: "<<right_coordination_exc<<" beside "<<sv_right_exc<<endl;

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

	//Biased_coordination
	if(vec_coordination.size())
	{
		++aside_influence;
		while(true)
		{
			uint vecSize(vec_coordination.size());
			if(vecSize == 1)
			{
				break;
			}

			map<uint, uint> map_coordination;
			for(register uint i(0); i < vecSize; ++i)
			{
				++map_coordination[vec_coordination[i].first + vec_coordination[i].second];
			}

			for(map<uint, uint>::iterator it(map_coordination.begin()); it != map_coordination.end(); ++it)
			{
				if((it->second * 100 / vecSize) > 50)
					gapped_read_count += it->second;
				else
					contradict_read_count += it->second;
			}

			break;
		}
	}


	//Low Complexity
	if(complexity_total / support_linkage_count < COMPLEXITY_THRESHOLD)
		sv_flag |= SV_LC;
	//GOSSIP("Marker 5 in [_ProcessSamCore]");
	//Aside influence
	if(aside_influence)
		sv_flag |= SV_AI;
	if(bias_influence)
		sv_flag |= SV_BIAS;

	//Miss linkage
	if(!support_linkage_count)
		sv_flag |= SV_ML;
	//GOSSIP("Marker 6 in [_ProcessSamCore]");
	//Support
	for(;;)
	{
		if((aligned_read_count + gapped_read_count) == 0)
		{
			sv_flag |= SV_F;
			break;
		}
		u16_t * sr(&gapped_read_count);
		u16_t * nsr(&aligned_read_count);
		if(*sr < SUPPORTING_READS_LOW)
			sv_flag |= SV_F;
		GOSSIP("Marker 7 in [_ProcessSamCore]");
		/*if(sv_type == STATE_I)
			NULL;
		else if(sv_type == STATE_D)
			swap(sr, nsr);
		else
			NR("NR point 2 in [_ProcessSamCore]");*/

		if(((*sr * 100 / (*sr + *nsr)) < 70) && ((*sr * 100 / (*sr + *nsr)) > 30))
			sv_flag |= SV_HET;

		if(*sr >= SUPPORTING_READS_LOW && *sr < SUPPORTING_READS_HIGH && ((*sr * 100 / (*sr + *nsr)) > 70) )
			sv_flag |= SV_LS;
		else if(*sr >= SUPPORTING_READS_HIGH && ((*sr * 100 / (*sr + *nsr)) > 70))
			sv_flag |= SV_HS;
		else
		{
			sv_flag |= SV_HET;
			sv_flag |= SV_F;
		}
		if(contradict_read_count * 100 / (contradict_read_count + aligned_read_count + gapped_read_count) > 50)
		{
			sv_flag &= (~SV_HS & ~SV_LS);
			sv_flag |= SV_F;
			sv_flag |= SV_BIAS;
		}
		//GOSSIP("Marker 8 in [_ProcessSamCore]");
		//cerr<<"[_ProcessSamCore]: Gapped/Aligned/Contradict: "<<gapped_read_count<<"/"<<aligned_read_count<<"/"<<contradict_read_count<<endl;
		break;
	}

	//OUTPUT
	string insertionGenotype;
	if(sv_type == STATE_I)
	{
		int cmpCount(0);
		for(unordered_map<string, uint>::iterator it(insertionCountingMap.begin()); it != insertionCountingMap.end(); ++it)
		{
			//cerr<<"[_ProcessSamCore]: Insertion alleles: "<<it->first<<"\t"<<it->second<<endl;
			if(it->second > cmpCount)
			{
				cmpCount = it->second;
				insertionGenotype = it->first;
			}
		}
		info = STATE_I;
		info += insertionGenotype;
	}
	if(((!(sv_flag & SV_HS)) && (!(sv_flag & SV_LS))))
		info.clear();
	//else if(sv_flag & SV_AI)
		//info.clear();


	info += "/";
	if(sv_flag & SV_F)
		info += "F";
	if(sv_flag & SV_LS)
		info += "L";
	if(sv_flag & SV_HS)
		info += "H";
	if(sv_flag & SV_ML)
		info += "M";
	if(sv_flag & SV_LC)
		info += "C";
	if(sv_flag & SV_AI)
		info += "A";
	if(sv_flag & SV_BIAS)
		info += "B";
	if(sv_flag & SV_HET)
	{
		info += "Z";
		if(sv_flag & SV_F)
		{
			info += sv_type;
			if(sv_type == STATE_D)
			{
				ushort tmpLen(p->len);
				info += lexical_cast<string>(tmpLen);
			}
			else if(sv_type == STATE_I)
				info += insertionGenotype;
		}
	}
	if(sv_type == STATE_D && ((sv_flag & SV_LS) || (sv_flag & SV_HS)))
	{
		info += "D(";
		info += sv_genotype;
		info += ")";
	}

	return (void*) NULL;
}

int ImportSAMOfWholeScaffold(MAP_UINT_STR & mapSam, string name)
{
	int count(0);
	//Import SAM result of a scaffold
	char isam_cmd[BUFF_LEN];
	name = SuitFilenameToShell(name);
	uint limit(bam_path.size());
	for(uint i(0); i < limit; ++i)
	{
		snprintf(isam_cmd, BUFF_LEN, "%s view %s %s 2>/dev/null", samtools_path, bam_path[i].c_str(), name.c_str());
		//cerr<<"ImportSAMOfWholeScaffold: "<<isam_cmd<<endl;
		DBG(isam_cmd);
		FILE * ISAM = popen(isam_cmd, "r");

		char buff[BUFF_LEN];

		string position;
		while(true)
		{
			fgets(buff, (BUFF_LEN - 1), ISAM);
			if(feof(ISAM))
				break;
			if(strlen(buff) == 0)
				continue;
			stringstream ss(buff);
			ss>>position>>position>>position>>position;
			GOSSIP(position);

			mapSam[lexical_cast<int>(position)] = buff;
			++count;
		}
		pclose(ISAM);
	}

	return count;
}

int SearchError(MAP_VEC_FASTA_BP & data, MAP_VEC_UINT & suspiciousLoci, int cpu)
{
	__SamProcessingSlot = __DetermineErrorsFromSam;
	ASSERT((__SamProcessingSlot == __DetermineErrorsFromSam), "__SamProcessingSlot should be assigned first.");

	int count(0);

	MAP_VEC_UINT_ITER itLimit(suspiciousLoci.end());
	for(MAP_VEC_UINT_ITER it(suspiciousLoci.begin()); it != itLimit; ++it)
	{
		cerr<<"Searching errors on "<<it->first<<"..."<<endl;

		MAP_UINT_STR mapSam;
		cerr<<ImportSAMOfWholeScaffold(mapSam, it->first)<<" records of "<<it->first<<" imported."<<endl;

		vector<SEARCH_ERROR_PARAM> param_vector(it->second.size(), SEARCH_ERROR_PARAM(data, mapSam));
		CThreadManager threadMgr(cpu);
		threadMgr.SetTimer(0.001);
		{
			uint limit(it->second.size());
			for(register uint i(0); i < limit; ++i)
			{
				param_vector[i].set(it->first, it->second[i]);
				threadMgr.AddThread(_SearchError, (void *)(&param_vector[i]));
			}
		}

		threadMgr.Run();

		uint innerCount(0);
		MAP_VEC_FASTA_BP_ITER it2(data.find(it->first));
		ASSERT(it2 != data.end(), "it2 != data.end() in [SearchError]");

		{
			uint limit(it2->second.size());
			for(register uint i(1); i < limit; ++i)
			{
				if(it2->second[i].app)
				{
					++innerCount;
					//cerr<<"[SearchError]: A candidate error @"<<i<<" allele: "<<it2->second[i].ch<<endl;
				}
			}
		}

		cerr<<innerCount<<" candidate errors are found in "<<it->first<<"."<<endl;
	}

	cerr<<"Candidates refreshing..."<<endl;
	for(MAP_VEC_FASTA_BP_ITER it(data.begin()); it != data.end(); ++it)
	{
		uint count(0);
		uint snpCount(0);
		uint limit(it->second.size());
		for(register uint i(1); i < limit; ++i)
		{
			if(it->second[i].app)
			{
				DBG(string("Cigar: ") + it->second[i].app->info);
				CigarMergeAndSelect(it->second[i].app->info);
				DBG(string("CigarAfter: ") + it->second[i].app->info);
				if(it->second[i].app->info.empty())
				{
					++count;
					it->second[i].freeApp();
				}
			}
			if(it->second[i].snp)
			{
				DBG(string("Snp: ") + it->second[i].snp->info);
				CureSNP(it->second[i].snp->info);
				DBG(string("SnpAfter: ") + it->second[i].snp->info);
				if(it->second[i].snp->info.empty())
				{
					it->second[i].freeSnp();
				}
				else
				{
					++snpCount;
					//Unmask next sentences to allow SNP modification inside the FASTA
					//it->second[i].ch = it->second[i].snp->info[0];
				}
			}
		}
		if(count)
			cerr<<count<<" errors eliminated due to low support."<<endl;
		if(snpCount)
			cerr<<snpCount<<" SNP curated."<<endl;
	}

	for(MAP_VEC_FASTA_BP_ITER it(data.begin()); it != data.end(); ++it)
	{
		uint limit(it->second.size());
		for(register uint i(1); i < limit; ++i)
		{
			if(it->second[i].app)
				++count;
		}
	}

	cerr<<"Totally "<<count<<" candidate errors are found."<<endl;

	return count;
}

void * _SearchError(void * param)
{
	SEARCH_ERROR_PARAM * p = (SEARCH_ERROR_PARAM*) param;

	string name(p->name);
	char buff[BUFF_LEN];

	snprintf(buff, BUFF_LEN, "%s%s_%u.fa", tmp_path, p->name.c_str(), p->locus);
	p->fasta_fn = strdup(buff);
	DBG(p->fasta_fn);
	//p->fasta_fn = (string(tmp_path) + p->name + lexical_cast<string>(p->locus) + ".fa");
	snprintf(buff, BUFF_LEN, "%s%s_%u.fq", tmp_path, p->name.c_str(), p->locus);
	p->fastq_fn = strdup(buff);
	DBG(p->fastq_fn);
	//p->fastq_fn = (string(tmp_path) + p->name + lexical_cast<string>(p->locus) + ".fq");
	snprintf(buff, BUFF_LEN, "%s%s_%u.sam", tmp_path, p->name.c_str(), p->locus);
	p->sam_fn = strdup(buff);
	DBG(p->sam_fn);
	//p->sam_fn = (string(tmp_path) + p->name + lexical_cast<string>(p->locus) + ".sam");
	string rm(p->fasta_fn);


	//Check remaining disk space
	struct statfs diskSpace;
	statfs(tmp_path, &diskSpace);
	if(diskSpace.f_bfree * diskSpace.f_bsize < 1024000)
	{
		cerr<<"No free space in temporaray folder! Program terminated."<<endl;
		exit(EXIT_FAILURE);
	}

	//cerr<<"__FastaPipe"<<endl;
	__FastaPipe(param);
	//cerr<<"__FastqPipe"<<endl;
	__FastqPipe(param);
	if(!p->fasta_bad_bit && !p->fastq_bad_bit)
	{
		//cerr<<"__SamPipe"<<endl;
		__SamPipe(param);
		ASSERT(__SamProcessingSlot, "__SamProcessingSlot should be defined at first in [_SearchError]")
		__SamProcessingSlot(param);
		remove(p->sam_fn);
	}
	else
	{
		DBG(p->name + ' ' + lexical_cast<string>(p->locus) + " bypassed.");
	}

	remove(p->fastq_fn);
	remove(p->fasta_fn);

#define RM(suffix) ((rm + suffix).c_str())
	remove(RM(".amb"));
	remove(RM(".ann"));
	remove(RM(".bwt"));
	remove(RM(".pac"));
	remove(RM(".rbwt"));
	remove(RM(".rpac"));
	remove(RM(".rsa"));
	remove(RM(".sa"));
#undef RM

	return (void*) NULL;
}

void * __FastaPipe(void * param)
{
	SEARCH_ERROR_PARAM * p = (SEARCH_ERROR_PARAM*) param;

	ofstream OFA(p->fasta_fn);
	if(!OFA)
	{
		cerr<<"Error opening "<<p->fasta_fn<<endl;
		p->fasta_bad_bit = 1;
		return (void*) NULL;
	}

	int left((int)p->locus);
	left -= interval;
	if(left < 0)
		left = 0;
//cerr<<"__FastaPipe: "<<left<<' '<<right<<endl;

	string fastaSeq;
	_ExportFastaSeq(p->data, p->name, (uint)left, (p->locus + p->len + interval), fastaSeq, true);

	if(fastaSeq.empty())
		p->fasta_bad_bit = 1;

	OFA<<">1\n"<<fastaSeq<<"\n\0"<<flush;
	OFA.close();

	return (void*) NULL;
}

void * __FastqPipe(void * param)
{
	SEARCH_ERROR_PARAM * p = (SEARCH_ERROR_PARAM*) param;

	ofstream OFQ(p->fastq_fn);
	if(!OFQ)
	{
		cerr<<"Error opening "<<p->fastq_fn<<endl;
		p->fastq_bad_bit = 1;
		return (void*) NULL;
	}

	uint count(0);
	int left(((int)p->locus + 1 - IS - 3*SD));
	if(left < 1)
		left = 1;
	uint right(p->locus + 2 + IS + 3*SD);

	//cerr<<"[_FastqPipe]: map_range: locus left right "<<p->locus<<' '<<left<<' '<<right<<endl;

	MAP_UINT_STR_ITER itHead(p->mapSam.lower_bound(left));
	MAP_UINT_STR_ITER itTail(p->mapSam.lower_bound(right));

	left = (p->locus - interval);
	right = (p->locus + interval);

	//cerr<<"[_FastqPipe]: capture_range: locus left right "<<p->locus<<' '<<left<<' '<<right<<endl;

	string tmpStr;
	for(; itHead != itTail; ++itHead)
	{
		tmpStr = (itHead->second);
		if(tmpStr.empty())
			continue;
		//cerr<<tmpStr<<endl;
		VSTR tmpVec;
		regex_split(back_inserter(tmpVec), tmpStr);

		uint position(lexical_cast<uint>(tmpVec[pos]));
		//cerr<<"[_FastqPipe]: Read @ position found: "<<position<<endl;
		if(tmpVec[cigar].find('*') == string::npos) //Read mapped
		{
			if((position < right) && ((position + tmpVec[seq].size()) >= left)) // Read is overlapped with target region
			{
				OFQ<<"@"<<count++<<"_"<<position<<"\n"<<tmpVec[seq]<<"\n+\n"<<tmpVec[qual]<<"\n";
				//cerr<<"[_FastqPipe]: Read @ position outputed 1: "<<position<<endl;
			}
		}
		else //Partly anchored
		{
			bool rev = (lexical_cast<u16_t>(tmpVec[flag]) & 0x20);
			uint temp1(0);
			uint temp2(0);
			if(!rev)   //F strand.
			{
				temp1 = position + IS - tmpVec[seq].size() - 3 * SD;
				temp2 = position + IS + 3 * SD;
			}
			else       //R strands
			{
				temp1 = position + tmpVec[seq].size() - IS - 3 * SD;
				temp2 = position + 2 * tmpVec[seq].size() - IS + 3 * SD;
			}
			if(temp1 < right && temp2 >= (uint)left)
			{
				OFQ<<"@"<<count++<<"_"<<position<<"\n"<<tmpVec[seq]<<"\n+\n"<<tmpVec[qual]<<"\n";
				//cerr<<"[_FastqPipe]: Read @ position outputed 2: "<<position<<endl;
			}
		}
	}
	OFQ<<"\0"<<flush;
	OFQ.close();

	if(!count)
		p->fastq_bad_bit = 1;
	//cerr<<count<<" in range read in [__FastqPipe]"<<endl;

	return (void*) NULL;
}

void * __SamPipe(void * param)
{
	SEARCH_ERROR_PARAM * p = (SEARCH_ERROR_PARAM*) param;

	char sys[BUFF_LEN];
	snprintf(sys, BUFF_LEN, "%s index %s 2>/dev/null; %s aln -o 1 -e 50 -L %s %s 2>/dev/null | %s samse %s - %s 1>%s 2>/dev/null",\
				bwa_path, SuitFilenameToShell(p->fasta_fn).c_str(), bwa_path, SuitFilenameToShell(p->fasta_fn).c_str(),\
				SuitFilenameToShell(p->fastq_fn).c_str(),  bwa_path, SuitFilenameToShell(p->fasta_fn).c_str(),\
				SuitFilenameToShell(p->fastq_fn).c_str(), SuitFilenameToShell(p->sam_fn).c_str());
	DBG(sys);
	if(system(sys))
		cerr<<"__SamPipe error on "<<p->name<<" "<<p->locus<<endl;

	return (void*) NULL;
}

void * __DetermineErrorsFromSam(void * param)
{
	SEARCH_ERROR_PARAM * p = (SEARCH_ERROR_PARAM*) param;

	ifstream I(p->sam_fn);
	ASSERT(I, "SAM input opening error!");
	string sam_content;
	MAP_VEC_FASTA_BP_ITER it(p->data.find(p->name));
	ASSERT((it != p->data.end()));
	uint zeroCoordination((p->locus - interval));
	uint leftValid_inc((p->locus - halfInterval));
	uint rightValid_exc((p->locus + halfInterval));

	//cerr<<"[__DetermineErrorsFromSam]: interval/halfinterval/locus/zeroCoordination: "<<interval<<'/'<<halfInterval<<'/'<<p->locus<<'/'<<zeroCoordination<<endl;

	while(true)
	{
		getline(I, sam_content, '\n');
		if(!I)
			break;
		VSTR line_vec;
		if(!regex_split(back_inserter(line_vec), sam_content))
		{
			DBG(string("Split Error: ") + sam_content);
			continue;
		}

		GOSSIP("Marker 1 at [__DetermineErrorsFromSam]");
		ASSERT((!line_vec[0].empty()), "The head column shouldn't be empty!");

		if(line_vec[0][0] == '@')
			continue;
		else if(!(atoi(line_vec[flag].c_str()) & SAM_UNMAPPED))
		{
			//u16_t rd_length(line_vec[seq].size());
			u32_t rd_pos(lexical_cast<u32_t>(line_vec[pos]) - 1);

			GOSSIP("Marker 2 at [__DetermineErrorsFromSam]");

			vector<pair<uint, char> > cigar_vec;
			__SeperateCigarPure(line_vec[cigar], cigar_vec);

			GOSSIP("Marker 3 at [__DetermineErrorsFromSam]");
			if(cigar_vec.size() != 1)
			{
				uint accumulation(0);
				uint limit(cigar_vec.size() - 1);
				for(register uint i(0); i < limit; ++i)
				{
					u32_t real_left_coordination_inc(zeroCoordination + rd_pos + accumulation);
					u32_t real_right_coordination_exc(real_left_coordination_inc + cigar_vec[i].first);
					if(cigar_vec[i].second == STATE_D)
					{
						GOSSIP("Marker 4 at [__DetermineErrorsFromSam]");
						if((real_left_coordination_inc >= leftValid_inc) && (real_right_coordination_exc < rightValid_exc) && (cigar_vec[i].first <= GAP_LIMIT))
						{
							//cerr<<"[__DetermineErrorsFromSam]: Candidate deletion "<<STATE_D<<lexical_cast<string>(cigar_vec[i].first)<<" @ "<<real_left_coordination_inc<<endl;
							it->second[real_left_coordination_inc].MutexOn();
							it->second[real_left_coordination_inc].InitApp()->info += STATE_D;
							it->second[real_left_coordination_inc].app->info += lexical_cast<string>(cigar_vec[i].first);
							it->second[real_left_coordination_inc].MutexOff();
						}
					}
					else if(cigar_vec[i].second == STATE_I)
					{
						GOSSIP("Marker 5 at [__DetermineErrorsFromSam]");
						accumulation += cigar_vec[i].first;
						if((real_left_coordination_inc >= leftValid_inc) && (real_right_coordination_exc < rightValid_exc) && (cigar_vec[i].first <= GAP_LIMIT))
						{
							//cerr<<"[__DetermineErrorsFromSam]: Candidate insertion "<<STATE_I<<lexical_cast<string>(cigar_vec[i].first)<<" @ "<<real_left_coordination_inc<<endl;
							it->second[real_left_coordination_inc].MutexOn();
							it->second[real_left_coordination_inc].InitApp()->info += STATE_I;
							it->second[real_left_coordination_inc].app->info += lexical_cast<string>(cigar_vec[i].first);
							it->second[real_left_coordination_inc].MutexOff();
						}
					}
					else
					{
						accumulation += cigar_vec[i].first;
						continue;
					}
				}
			}

			if(snpCorrection)
			{
				ASSERT(line_vec.size()>=12, "Column should larger than 12 when aligned in [__DetermineErrorsFromSam]");

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
						ASSERT(line_vec[j].size()>2, "Tags should have more than 2 characters in [__DetermineErrorsFromSam]");

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
							u32_t snp_coordination(zeroCoordination + rd_pos + clipAddon + cigarVec[i].first);
							GOSSIP("Marker 6 at [__DetermineErrorsFromSam]");
							if((snp_coordination >= leftValid_inc) && (snp_coordination < rightValid_exc))
							{
								it->second[snp_coordination].MutexOn();
								it->second[snp_coordination].InitSnp()->info += cigarVec[i].second;
								it->second[snp_coordination].MutexOff();
							}
						}
					}
				}
			}
		}
	}

	return (void *) NULL;
}

inline void CigarMergeAndSelect(string & str)
{
	ASSERT(str.size(), "Input cigar string shouldn't be empty in [CigarMergeAndSelect]");
	typedef unordered_map<string, uint> MAP_STR_UINT;
	MAP_STR_UINT cigar_map;
	string cigar;
	cigar += str[0];
	uint limit(str.size());
	for(register uint i(1); i < limit; ++i)
	{
		if(isalpha(str[i]))
		{
			++cigar_map[cigar];
			cigar.clear();
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
		if(it->second > count)
		{
			str = it->first;
			count = it->second;
		}
		else if(it->second == count)
		{
			if(str[0] == STATE_I && it->first[0] == STATE_D)
				str = it->first;
		}
	}
	//cerr<<"[CigarMergeAndSelect]: FinalCigar: "<<str<<endl;
}

inline void CureSNP(string & str)
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

int Base2BaseCuration(MAP_VEC_FASTA_BP & data, int cpu)
{
	__SamProcessingSlot = _ProcessSamCore;
	ASSERT((__SamProcessingSlot == _ProcessSamCore), "__SamProcessingSlot should be assigned first.");

	int count(0);

	for(MAP_VEC_FASTA_BP_ITER it(data.begin()); it != data.end(); ++it)
	{
		if(it->second.size() < uint(2 * interval))
		{
			cerr<<it->first<<" too short, bypassed."<<endl;
			continue;
		}

		cerr<<"Performing Base2base curation on "<<it->first<<endl;

		uint limit(it->second.size());
		CThreadManager threadMgr(cpu);
		threadMgr.SetTimer(0.001);

		uint innerCount(0);
		for(register uint i(1); i < limit; ++i)
		{
			if(it->second[i].app)
			{
				++innerCount;
			}
		}

		MAP_UINT_STR mapSam;
		cerr<<ImportSAMOfWholeScaffold(mapSam, it->first)<<" records of "<<it->first<<" imported."<<endl;

		vector<SEARCH_ERROR_PARAM> param_vector(innerCount, SEARCH_ERROR_PARAM(data, mapSam));
		count += innerCount;
		cerr<<innerCount<<" candidate errors are exist in "<<it->first<<"."<<endl;
		cerr<<"Processing..."<<endl;
		string lenStr;
		for(register uint i(1),j(0); i < limit && j < innerCount; ++i)
		{
			if(it->second[i].app)
			{
				param_vector[j].set(it->first, i);
				ASSERT(!it->second[i].app->info.empty(), "it->second[i].app->info.empty() in [Base2BaseCuration]");
				try {
					param_vector[j].len = lexical_cast<ushort>(it->second[i].app->info.substr(1, string::npos));
				}
				catch (exception_detail::clone_impl<exception_detail::error_info_injector<bad_lexical_cast> > & e) {
					cerr<<"Bad lexical cast in Base2BaseCuration: "<<endl
						<<e.what()<<endl
						<<"it->second[i].app->info: "<<it->second[i].app->info<<endl;
				}
				catch (out_of_range & e) {
					cerr<<"Out of range in Base2BaseCuration: "<<endl
						<<e.what()<<endl
						<<"it->second[i].app->info: "<<it->second[i].app->info<<endl;
				}
				ASSERT(param_vector[j].len, "param_vector[j++].len, length shouldn't be 0 in [Base2BaseCuration]");
				threadMgr.AddThread(_SearchError, (void*)&param_vector[j]);
				++j;
			}
		}
		threadMgr.Run();

		uint bypassedCount(0);
		MAP_VEC_FASTA_BP_ITER it2(data.find(it->first));
		ASSERT(it2 != data.end(), "it2 != data.end() in [SearchError]");

		{
			uint previousMark(0);
			uint accumulatedError(0);

			uint limit(it2->second.size());
			for(register uint i(1); i < limit; ++i)
			{
				if(it2->second[i].app)
				{
					if((i - previousMark <= DETERMINE_CONFUSION_STEP) && (i > DETERMINE_CONFUSION_STEP))
					{
						++accumulatedError;
					}
					else
					{
						if(accumulatedError)
						{
							for(register uint j(previousMark); j < i; ++j)
							{
								if(it2->second[j].app)
								{
									//cerr<<"[Base2BaseCuration]: Correction @ "<<j<<'\t'<<it->second[j].app->info<<" eliminated"<<endl;
									it2->second[j].freeApp();
								}
							}
							bypassedCount += accumulatedError;
						}
						accumulatedError = 0;
						previousMark = i;
					}
				}
			}
		}

		if(bypassedCount)
			cerr<<bypassedCount<<" corrected result(s) eliminated in "<<it->first<<" due to confusion."<<endl;

		cerr<<it->first<<" curated."<<endl;
	}

	return count;
}

int EliminateError(MAP_VEC_FASTA_BP & data)
{
	int count(0);

	for(MAP_VEC_FASTA_BP_ITER it(data.begin()); it != data.end(); ++it)
	{
		uint limit(it->second.size());
		for(uint i(0); i < limit; ++i)
		{
			if(it->second[i].app)
			{
				ASSERT((it->second[i].app->info.size() >= 2), "it->second[i].app->info.size() shouldn't be <2 at [Base2BaseCuration]");
				if(it->second[i].app->info[0] == STATE_I)
				{
					if(isdigit(it->second[i].app->info[1]))
					{
						++count;
						it->second[i].freeApp();
					}
				}
			}
		}
	}

	return count;
}

//Returns the length of each type in vector
uint __SeperateCigarPure(string & str, VEC_UINT_CHAR & vec)
{
	string pos;
	uint limit(str.size());
	for(register uint i(0); i < limit; ++i)
	{
		if(isalpha(str[i]))
		{
			ASSERT(pos.size());
			vec.push_back(make_pair((atoi(pos.c_str())), str[i]));
			pos.clear();
		}
		else
			pos += str[i];
	}

	return vec.size();
}

uint __SeperateCigarMismatch(string & str, VEC_UINT_CHAR & vec)
{
	string pos;
	uint character(0);
	uint accumulate(0);
	uint limit(str.size());
	for(register uint i(0); i < limit; ++i)
	{
		if(isalpha(str[i]))
		{
			if(!pos.empty())
				accumulate += atoi(pos.c_str());
			vec.push_back(make_pair(accumulate + character, str[i]));
			pos.clear();
			++character;
		}
		else if(str[i] == '^')
		{
			if(!pos.empty())
				accumulate += atoi(pos.c_str());
			pos.clear();
			++i;
			while(isalpha(str[i]))
			{
				++i;
				if(i == limit)
					break;
			}
			--i;
		}
		else
			pos += str[i];
	}

	return character;
}

int SimulateIndel(MAP_VEC_FASTA_BP & data)
{
	int count(0);

	vector<uint> indelLengthVec;
	for(uint i(INDEL_LOWER_LIMIT); i < INDEL_UPPER_LIMIT; ++i)
	{
		uint limit(10000.0 / ((double)i * (double)i));
		for(register uint j(0); j < limit; ++j)
		{
			indelLengthVec.push_back(i);
		}
	}
	double indelLengthVecSize(indelLengthVec.size());
	cerr<<"Indel length Vector size: "<<indelLengthVecSize<<endl;

	vector<uint> indelQuantityForEach;
	uint totalLength(0);
	for(MAP_VEC_FASTA_BP_ITER it(data.begin()); it != data.end(); ++it)
	{
		indelQuantityForEach.push_back(it->second.size());
		totalLength += it->second.size();
	}
	cerr<<"Total Length of all sequences: "<<totalLength<<endl;

	if(INDEL_QUANTITY > (totalLength * 0.1))
	{
		cerr<<"Defined indel quality has been set to (Total length of sequence * 0.1)"<<endl;
		INDEL_QUANTITY = totalLength * 0.1;
	}

	for(uint i(0); i < indelQuantityForEach.size(); ++i)
	{
		//cerr<<indelQuantityForEach[i]<<'\t'<<INDEL_QUANTITY<<'\t'<<totalLength<<'\t'<< indelQuantityForEach[i] * INDEL_QUANTITY / totalLength<<endl;
		indelQuantityForEach[i] = (double)indelQuantityForEach[i] * (double)INDEL_QUANTITY / (double)totalLength;
	}

	int i(0);
	for(MAP_VEC_FASTA_BP_ITER it(data.begin()); it != data.end(); ++it)
	{
		srand48(time(0));
		uint limit(indelQuantityForEach[i++]);
		cerr<<"Indel quantity for "<<it->first<<": "<<limit<<endl;
		for(uint j(0); j < limit; )
		{
			uint locus((int)((double)it->second.size() * (drand48())));
			if(!it->second[locus].app)
			{
				++count;
				++j;
				it->second[locus].InitApp();
				string & str(it->second[locus].app->info);
				char ch("ID"[(int)(2.0 * drand48())]);
				int length(indelLengthVec[(int)(indelLengthVecSize * (drand48()))]);
				if(ch == STATE_D)
				{
					str += ch;
					str += lexical_cast<string>(length);
				}
				else if(ch == STATE_I)
				{
					str += ch;
					for(register int i(0); i < length; ++i)
					{
						str += "ATGC"[(int)(4.0 * drand48())];
					}
				}
				//cerr<<"[SimulateIndel]: Indel "<<str<<" at "<<locus<<endl;
			}
		}
	}

	if(HETEROZYGOUS_RATE)
	{
		srand48(time(0));
		cerr<<"[SimulateIndel]: Heterozygous mode on, rate: "<<HETEROZYGOUS_RATE<<endl;
		vector<MAP_VEC_FASTA_BP_ITER> iterVec4AllSeq;
		for(MAP_VEC_FASTA_BP_ITER it(data.begin()); it != data.end(); ++it)
		{
			iterVec4AllSeq.push_back(it);
		}
		uint totalHetSite(ceil(totalLength * HETEROZYGOUS_RATE / 3 * 4));
		cerr<<"[SimulateIndel]: TotalHetSites: "<<totalHetSite<<endl;
		for( ; totalHetSite != 0; --totalHetSite)
		{
			VEC_FASTA_BP & currentVec(iterVec4AllSeq[iterVec4AllSeq.size() * drand48()]->second);
			uint locus(currentVec.size() * drand48());
			char orgChar(currentVec[locus].ch);
			currentVec[locus].ch = "ATGC"[(int)(4.0 * drand48())];
			//if(orgChar != currentVec[locus].ch)
			//	cerr<<"[SimulateIndel]: "<<orgChar<<"-> "<<currentVec[locus].ch<<" @ "<<locus<<endl;
		}
	}

	return count;
}

