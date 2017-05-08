# Here is the example of the programe commands of structures variation detection.

# BWA alignment
bwa bwasw -t 4 $HG19 $CONTIG


# SV detection by SOAPSV. We use SOAPSV 1.02 to call SVs. As the pipeline is sophisticated and the steps are trivial, we don't want to show all of them in this script, for more detials, please refer into the documents in the SOAPSV directory.


lastz --targetcapsule=$CAPSULE $FASTA[nameparse=darkspace] --strand=both --chain --ambiguous=iupac --gapped --ydrop=50000 --gap=1000,1 --format=axt --output=$AXT --markend
axtSort $AXT > $SORT_AXT
...
best_hit $SORT_AXT > $BEST_AXT
...
intro_indel_1.3 $FINAL_AXT > $SV

# SV detection by Pindel  
pindel -f   $HG19 -i $CFG -o $OUT_PREFIX 
pindel2vcf -P $PREFIX -r $HG19 -R hg19 -d hg19 -v $VCF

# SV detection by CNVnator
./cnvnator -genome hg19 -root out.root  -tree $BAM
./cnvnator -genome hg19 -root out.root -his 100
./cnvnator -root out.root  -stat 100
./cnvnator -root out.root -partition 100
./cnvnator -root out.root  -call 1000

# SV detection by breakdancer
bam2cfg.pl -q 20 -c 3 -g -h $BAM > $CFG
breakdancer -o $PREFIX -q 20 -d $CTX -a -y 30 $CFG

# SV detection by GenomeSTRiP(svtoolkit)
java -cp ${classpath} ${mx}     org.broadinstitute.sting.queue.QCommandLine     -S ${SV_DIR}/qscript/SVPreprocess.q     -S ${SV_DIR}/qscript/SVQScript.q     -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar     -cp ${classpath}     -configFile conf/genstrip_parameters.txt     -disableGATKTraversal     -tempDir ${SV_TMPDIR}     -R $HG19     -computeGCProfiles     -genomeMaskFile hg19.mask.101.fasta     -ploidyMapFile hg19.ploidy.map     -copyNumberMaskFile cn2_mask_hg19.fasta     -genderMapFile gender.list     -runDirectory ${runDir}     -computeGCProfiles     -md ${runDir}/metadata     -jobLogDir ${runDir}/logs     -I ${bam}     --disableJobReport     -run     || exit 1

java -cp ${classpath} ${mx}     org.broadinstitute.sting.queue.QCommandLine     -S ${SV_DIR}/qscript/SVDiscovery.q     -S ${SV_DIR}/qscript/SVQScript.q     -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar     --disableJobReport     -cp ${classpath}     -configFile ./genstrip_parameters.txt     -tempDir ${SV_TMPDIR}     -R $HG19     -genomeMaskFile hg19.mask.101.fasta     -genderMapFile gender.list     -runDirectory ${runDir}     -md ${runDir}/metadata     -disableGATKTraversal     -jobLogDir ${runDir}/logs     -minimumSize 50     -maximumSize 1000000     -windowSize 20000000     -windowPadding 10000     -I ${bam}     -O ${sites}     -P select.validateReadPairs:false     -run     || exit 1


# Follwing, we combine the deletions in individual level between different methods and in popluation level among different individuals with using in-house scrtips. The methods are similar in 1000 genome paper. Merge exact breakpoint by locations and merge imprecise breakpoint by confident region.

# Genotyping deletions by GenomeSTRiP(svtoolkit)

java  -cp ${classpath} ${mx}     org.broadinstitute.sting.queue.QCommandLine     -S ${SV_DIR}/qscript/SVGenotyper.q     -S ${SV_DIR}/qscript/SVQScript.q     -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar     --disableJobReport     -cp ${classpath}     -configFile genstrip_parameters.txt     -tempDir ${SV_TMPDIR}     -R $HG19     -genomeMaskFile hg19.mask.101.fasta     -genderMapFile gender.list     -runDirectory ${runDir}     -md ${runDir}/metadata     -jobLogDir ${runDir}/logs     -I ${bam}     -vcf ${sites}     -disableGATKTraversal     -O ${genotypes}     -run     || exit 1
