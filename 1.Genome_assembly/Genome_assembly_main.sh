# Here is the example of the programe commands of raw data preprocess and genome assembly.

# Raw data preprocess
# filter_data_gz and duplication are in-house programes, but we can use a updated and open-source version programe SOAP nuke, the website is https://github.com/BGI-flexlab/SOAPnuke
filter_data_gz -y -z -w 10 -B 40 -l 10827 -a 2 -b 3 -c 2 -d 3 1.fq.gz 2.fq.gz 1.fq.gz.reads.stat 1.fq.gz.clean.fq.gz  2.fq.gz.clean.fq.gz
duplication  1.fq.gz.clean 2.fq.gz.clean 1.fq.gz.clean.dup.clean 2.fq.gz.clean.dup.clean 1.fq.gz.clean.dup.stat
# kmer_freq_pfile and correct_error_pfile are in-house programes, but we can use a updated and open-source version programe ErrorCorrection, the website is http://sourceforge.net/projects/soapdenovo2/files/ErrorCorrection/SOAPec_v2.01.tar.gz
kmer_freq_pfile -k 17 -t 16 -p filterDupFq.list
correct_error_pfile -k 17 -l 6  -a 0 -t 16 -j 1 HG00684.freq.cz filterDupFq.list

# Soap Denovo assembly
# SOAPdenovo-63mer-V1.06_0718 is in-house programes, but we can use a updated and open-source version programe SOAPdenovo2, the website is https://github.com/aquaskyline/SOAPdenovo2
SOAPdenovo-63mer-V1.06_0718 pregraph -s $CFG -K 63 -o $PREFIX -p 16 
SOAPdenovo-63mer-V1.06_0718 contig -g $PREFIX -M 2 
SOAPdenovo-63mer-V1.06_0718 map -s $CFG -k 63 -g $PREFIX -p 16 
SOAPdenovo-63mer-V1.06_0718 scaff -g $PREFIX -p 16
