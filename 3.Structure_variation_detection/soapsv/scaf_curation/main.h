/*
 *  Project: SV judgment by unmapped mate alignment
 * Filename: main.h
 *
 *  Created on: Feb 8, 2010
 *      Author: Ruibang Luo
 *
 *     History:
 *         1.
 */

#pragma once
#ifndef MAIN_H_AQUA_
#define MAIN_H_AQUA_

#include<vector>
#include<string>

using namespace std;

void usage();
void usage(string);

//extern char * bam_path;
extern vector<string> bam_path;
extern char * samtools_path;
extern char * bwa_path;
extern char * tmp_path;
extern unsigned int readLength;
extern unsigned int IS;
extern unsigned int SD;
extern unsigned int INDEL_UPPER_LIMIT;
extern unsigned int INDEL_LOWER_LIMIT;
extern unsigned int INDEL_QUANTITY;
extern bool snpCorrection;
extern unsigned int SNP_LOWER_LIMIT;
extern unsigned int SUPPORTING_READS_LOW;
extern unsigned int SUPPORTING_READS_HIGH;
extern double HETEROZYGOUS_RATE;

#endif /* MAIN_H_ */
