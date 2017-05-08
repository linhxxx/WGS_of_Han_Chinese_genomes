/*
 *  Project: SV judgment by unmapped mate alignment
 * Filename: main.cpp
 *
 *  Created on: Feb 8, 2010
 *      Author: Ruibang Luo
 *
 *     History:
 *         1.
 */

#define _EXAM_ASSERT_TEST_
#include"main.h"
#include"general.h"
#include"func.h"
#include"func2.h"
#include"pileup.h"
#include<string>
#include<iostream>
#include<getopt.h>

using namespace std;
using namespace boost;
using namespace aqua;

//char * bam_path;
vector<string> bam_path;
char * samtools_path;
char * bwa_path;
char * tmp_path;
unsigned int readLength;
unsigned int IS;
unsigned int SD;
unsigned int INDEL_UPPER_LIMIT;
unsigned int INDEL_LOWER_LIMIT;
unsigned int INDEL_QUANTITY;
bool snpCorrection;
unsigned int SNP_LOWER_LIMIT;
unsigned int SUPPORTING_READS_LOW;
unsigned int SUPPORTING_READS_HIGH;
double HETEROZYGOUS_RATE;

int main(int argc, char ** argv)
{
	if(argc <= 1)
		usage();

	int cpu(4);
	int refBasedMode(0), scaffBasedMode(0), simulationMode(0), pileUpMode(0);
	snpCorrection = false;
	int globalCorrection(0);
	readLength = 0;
	char * input_fn(strdup("/dev/stdin"));
	char * output_fn(strdup("/dev/stdout"));
	char * output_fasta(NULL);
	char * suspicious_list(NULL);
	//bam_path = 0;
	bam_path.clear();
	samtools_path = strdup("/share/backup/luoruibang/bin/samtools");
	bwa_path = strdup("/share/backup/luoruibang/bin/bwa");
	tmp_path = strdup("/tmp/");
	IS = 200;
	SD = 10;
	INDEL_LOWER_LIMIT = 1; INDEL_UPPER_LIMIT = 500; INDEL_QUANTITY = 0;
	SNP_LOWER_LIMIT = 5;
	SUPPORTING_READS_LOW = 3; SUPPORTING_READS_HIGH = 10;
	HETEROZYGOUS_RATE = 0;
	int c;
	while((c = getopt(argc, argv, "H:p:i:o:hrsnP:kal:d:f:b:w:m:t:u:v:q:z:e:y:j:x:")) != -1)
	{
		//Already used parameters: "abcdefhijklmnopqrstuvwxyz HP"
		switch(c)
		{
		case 'r': refBasedMode = 1; break;
		case 's': scaffBasedMode = 1; break;
		case 'n': simulationMode = 1; break;
		case 'P': pileUpMode = atoi(optarg); break;
		case 'k': snpCorrection = true; break;
		case 'p': cpu = atoi(optarg); break;
		case 'i': free(input_fn); input_fn = strdup(optarg); break;
		case 'o': free(output_fn); output_fn = strdup(optarg); break;
		case 'f': output_fasta = strdup(optarg); break;
		case 'b':
			{
				--optind;
				do
				{
					//cerr<<argv[optind]<<" pushed."<<endl;
					bam_path.push_back(argv[optind++]);
				} while(optind < argc && strlen(argv[optind]) >= 1 && argv[optind][0] != '-');
				cerr<<bam_path.size()<<" sorted bam file path received."<<endl;
				break;
			}
		case 'a': globalCorrection = 1; break;
		case 'l': suspicious_list = strdup(optarg); break;
		case 'd': readLength = atoi(optarg); break;
		case 'w': free(bwa_path); bwa_path = strdup(optarg); break;
		case 'm': free(samtools_path); samtools_path = strdup(optarg); break;
		case 't': free(tmp_path); tmp_path = strdup(optarg); strcat(tmp_path, "/"); break;
		case 'u': IS = atoi(optarg); break;
		case 'v': SD = atoi(optarg); break;
		case 'H': HETEROZYGOUS_RATE = atof(optarg); break;
		case 'q': INDEL_LOWER_LIMIT = atoi(optarg); break;
		case 'z': INDEL_UPPER_LIMIT = atoi(optarg); break;
		case 'e': INDEL_QUANTITY = atoi(optarg); break;
		case 'y': SNP_LOWER_LIMIT = atoi(optarg); break;
		case 'j': SUPPORTING_READS_LOW = atoi(optarg); break;
		case 'x': SUPPORTING_READS_HIGH = atoi(optarg); break;
		case 'h': usage(); break;
		default: usage(argv[optind]); break;
		}
	}
	{
		if(SUPPORTING_READS_HIGH <= SUPPORTING_READS_LOW)
			usage("High threshold should larger than low threshold of supporing reads.");
		if(!(refBasedMode ^ scaffBasedMode ^ simulationMode ^ pileUpMode))
			usage("At least and only one mode should be selected.");
		if(scaffBasedMode)
		{
			if(!readLength)
				usage("Please set the shortest read length in libraries using -d.");
			if((globalCorrection && suspicious_list) || (!globalCorrection && !suspicious_list))
				usage("One of using global Correction or suspicious_list should be selected.");
			if(!input_fn)
				usage("Please use -i to import the assembly in fasta format.");
			if(!output_fasta)
				usage("Please use -f to export the assembly in fasta format.");
			if(/*!bam_path*/ bam_path.empty())
				usage("Please use -b to import the Merged BAM.");
			if(!readLength)
				usage("Read length not defined.");
			/*if((!FileExist(bam_path, false)))
				usage("Merged alignment file missing!");
			if((!FileExist((string(bam_path) + ".bai").c_str(), false)))
				usage("Merged alignment index missing, please use \"samtools index <merged bam>\" to create!");*/
			for(uint i(0); i < bam_path.size(); ++i)
			{
				if((!FileExist(bam_path[i].c_str(), false)))
					usage("Merged alignment file missing!");
				if((!FileExist((string(bam_path[i]) + ".bai").c_str(), false)))
					usage("Merged alignment index missing, please use \"samtools index <merged bam>\" to create!");
			}
			if((!FileExist(bwa_path, false)))
				usage("BWA missing, please use -w to specify a path.");
			if((!FileExist(samtools_path, false)))
				usage("SAMTOOLS missing, please use -w to specify a path.");
			if(!system((string("touch ") + tmp_path + "1").c_str()))
				system((string("rm ") + tmp_path + "1").c_str());
			else
			{
				free(tmp_path); tmp_path = strdup("/tmp/");
				if(!system((string("touch ") + tmp_path + "1").c_str()))
				{
					cerr<<"User defined temporary directory \""<<tmp_path<<"\" not available, using /tmp instead."<<endl;
					system((string("rm ") + tmp_path + "1").c_str());
				}
				else
					usage("Temporatory directory not available, please specify another one using -t.");
			}
		}
		if(simulationMode)
		{
			if(!input_fn)
				usage("Please use -i to import the sequences in fasta format.");
			if(!output_fasta)
				usage("Please use -f to export the sequences in fasta format.");
			if(!INDEL_QUANTITY && !HETEROZYGOUS_RATE)
				usage("Please use -e to define how many indels to generate.");
			if(INDEL_UPPER_LIMIT <= INDEL_LOWER_LIMIT)
				usage("INDEL_LOWER_LIMIT shouldn't be >= INDEL_UPPER_LIMIT.");
		}
		if(pileUpMode)
		{
			if(pileUpMode > (PILEUP_MODE_INDEL | PILEUP_MODE_SNP) || pileUpMode < 0x1)
			{
				usage("Pileup mode error.");
			}
			for(uint i(0); i < bam_path.size(); ++i)
			{
				if((!FileExist(bam_path[i].c_str(), false)))
					usage("Merged alignment file missing!");
			}
		}
	}


	if(refBasedMode)
	{
		cerr<<"Reference based mode selected."<<endl;
		VSTR file_list_vec;
		cerr<<"Importing sam list..."<<endl;
		ImportSamList(input_fn, file_list_vec);
		cerr<<"Processing sam list..."<<endl;
		ProcessSamList(file_list_vec, cpu, output_fn);
	}
	else if(scaffBasedMode)
	{
		cerr<<endl<<"Scaffold based mode selected."<<endl;
		MAP_VEC_FASTA_BP data;
		cerr<<endl<<"Importing Fasta..."<<endl;
		ImportFasta(input_fn, data, true);

		MAP_VEC_UINT suspiciousLoci;
		if(suspicious_list)
		{
			cerr<<endl<<"Importing suspicious list..."<<endl;
			ImportSuspiciousList(suspicious_list, data, true, readLength, suspiciousLoci);
		}
		else
		{
			cerr<<endl<<"Globally determining interval in sequences..."<<endl;
			GlobalDetermineSuspicious(data,readLength, suspiciousLoci);
		}
		cerr<<endl<<"Locating suspicious loci in sequences..."<<endl;
		SearchError(data, suspiciousLoci, cpu);
		cerr<<endl<<"Performing Base2Base curation..."<<endl;
		Base2BaseCuration(data, cpu);
		EliminateError(data);
		cerr<<endl<<"Exporting Fasta..."<<endl;
		ExportFasta(output_fasta, data, true);
		cerr<<endl<<"Exporting Details..." <<endl;
		ExportDetails(output_fn, data, true);
	}
	else if(simulationMode)
	{
		cerr<<endl<<"Simulation mode selected."<<endl;
		MAP_VEC_FASTA_BP data;
		cerr<<endl<<"Importing Fasta..."<<endl;
		ImportFasta(input_fn, data, true);
		cerr<<endl<<"Simulating indels..."<<endl;
		SimulateIndel(data);
		cerr<<endl<<"Exporting Fasta..."<<endl;
		ExportFasta(output_fasta, data, true);
		cerr<<endl<<"Exporting Details..." <<endl;
		ExportDetails(output_fn, data, true);
	}
	else if(pileUpMode)
	{
		int count(0);
		cerr<<"Pileup mode seleted"<<endl;
		MSMUS indel, snp;
		cerr<<"Starting pileup..."<<endl;
		count = Pileup(pileUpMode, indel, snp);
		cerr<<count<<" indel reports collected."<<endl;
		cerr<<"Polishing results..."<<endl;
		count = PileupPolish(pileUpMode, indel, snp);
		cerr<<count<<" sites polished"<<endl;
		cerr<<"Outputing results..."<<endl;
		count = PileupOutput(pileUpMode, indel, snp, output_fn);
		cerr<<count<<" indel results outputed."<<endl;
	}

	cerr<<endl<<"Job finished."<<endl;

	return 0;
}

void usage()
{
	cerr<<line
		<<"Assembly Restrictor"<<"\n"
		<<"Version 1.3"<<"\n"
		<<"Author: Ruibang Luo"<<"\n"
		<<"\n"
		<<"Parameters:"<<"\n"
		<<"-p [int]        Threads"<<"\n"
		<<"\n"
		<<"Reference based:"<<"\n"
		<<"-r              Enable this mode"<<"\n"
		<<"-i [string]     Input SAM files list"<<"\n"
		<<"-o [string]     Output result file"<<"\n"
		<<"\n"
		<<"Assembly based:"<<"\n"
		<<"-s              Enable this mode"<<"\n"
		<<"-a              Global Correction"<<"\n"
		<<"-k              Enable SNP correction"<<"\n"
		<<"-l [string]     Suspicious list, 1 as leftmost. Format: [seq_id start(include) end(exclude) | seq_id locus]"<<"\n"
		<<"-d [int]        The Length of the shortest read in all libraries"<<"\n"
		<<"-i [string]     Input the scafSeq"<<"\n"
		<<"-o [string]     Output result file, Default: /dev/stdout"<<"\n"
		<<"-f [string]     Output corrected fasta"<<"\n"
		<<"-b [strings]    Merged sorted BAM"<<"\n"
		<<"-w [string]     Specify the BWA program path, Default: /share/backup/luoruibang/bin/samtools"<<"\n"
		<<"-m [string]     Specify the SAMTOOLS program path, Default: /share/backup/luoruibang/bin/bwa"<<"\n"
		<<"-t [string]     Specify the temperary directory, Default: /tmp/"<<"\n"
		<<"-u [int]        Specify the insert size, Default: 200"<<"\n"
		<<"-v [int]        Specify the standard division, Default: 10"<<"\n"
		<<"-y [int]        Specify the lower limit of reads support for a SNP alternation, Default: 5"<<"\n"
		<<"-j [int]        Specify minimum supporting read(s) for a correction. Default: 3"<<"\n"
		<<"-x [int]        Specify how many supporting read(s) to define a high confident correction. Default: 10"<<"\n"
		<<"\n"
		<<"Simulation:"<<"\n"
		<<"-n              Enable this mode"<<"\n"
		<<"-i [string]     Input the scafSeq"<<"\n"
		<<"-o [string]     Output result file, Default: /dev/stdout"<<"\n"
		<<"-f [string]     Output corrected fasta"<<"\n"
		<<"-q [int]        Set the lower limitation of length of simulated indels, Default: 1"<<"\n"
		<<"-z [int]        Set the upper limitation of length of simulated indels, Default: 500"<<"\n"
		<<"-e [int]        How many indels to generate"<<"\n"
		<<"-H [double]     Simulate heterozygous rate of the genome"<<"\n"
		<<"\n"
		<<"Pileup:"<<"\n"
		<<"-P [int]        Enable this mode, 1 for InDel, 2 for SNP, 3 for both"<<"\n"
		<<"-o [string]     Output result file, Default: /dev/stdout"<<"\n"
		<<"-b [strings]    Merged sorted BAM"<<"\n"
		<<"-y [int]        Specify the lower limit of reads support for a SNP alternation, Default: 5"<<"\n"
		<<"-j [int]        Specify minimum supporting read(s) for a correction. Default: 3"<<"\n"
		<<"-h              This help"
		<<line;

	exit(EXIT_FAILURE);
}

void usage(string str)
{
	cerr<<str<<endl;
	usage();
}
