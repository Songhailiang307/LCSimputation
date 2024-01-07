#!/bin/bash
#PBS -q unlimited
#PBS -l nodes=1:ppn=2
#PBS -N gatk+beagle5

## "PBS" routine
if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
    echo [`date`]INFO: qsub routine env!
    echo [`date`]WORKDIR: $PBS_O_WORKDIR
    ncpus=`wc -l < $PBS_NODEFILE`
else
    echo [`date`]INFO: shell routine env!
    echo [`date`]INFO: working space [`pwd`]
    filename=$1
    chr=$2
fi

# set the header of fastq files
if [ ! -n "$filename" ]; then
    echo [`date`]ERROR: please specify the header of the fastq file!
    echo [`date`]EXAMPLE: 
    echo '    pair_end_1: /home/jiangyf/test/TEST/hello/test/S483403_H3JGGDMXX_L1_1.clean.fq.gz' 
    echo '    pair_end_2: /home/jiangyf/test/TEST/hello/test/S483403_H3JGGDMXX_L1_2.clean.fq.gz'
    echo '           OUT: S483403'
    echo '           DIR: /home/jiangyf/test/TEST/hello/test'
    echo '    qsub -v filename=/home/jiangyf/test/TEST/hello/test/S483403_H3JGGDMXX_L1_1.clean.fq.gz ./fastq2bam.PBS.sra'
    exit 1
else
    echo [`date`]INFO: file name [$filename]
    DIR=`dirname $filename`
    FILE=`basename $filename .list`
    cd $DIR
    echo [`date`]INFO: working space [`pwd`]
fi


# set files path that is needed in analysis
REF_FASTA_FILE=/home/songhailiang/Novogene_analysis/reference_download/Aru_download.fa  #this is the reference genome from web.
TEMP_DIR=/data/songhl/temp/
echo [`date`]REFERENCE: [$REF_FASTA_FILE]

# set software path
trim_path=/home/songhailiang/software/quality-trim
JAVA=java
GATK=/home/songhailiang/software/GATK/GenomeAnalysisTK.jar
PICARD=/home/songhailiang/software/GATK/picard-2.20.6/picard.jar
IlluQC=/home/songhailiang/software/NGSQCToolkit_v2.3.3/QC/IlluQC.pl
Trimming=/home/songhailiang/software/NGSQCToolkit_v2.3.3/Trimming/TrimmingReads.pl
FASTQC=/home/songhailiang/software/fastqc/FastQC/fastqc
BWA=/home/songhailiang/software/bwa/bwa-0.7.17/bwa
SAMTOOLS=/home/songhailiang/software/SAM_tools/samtools-1.9/samtools
# set routine parameters
max_mem=40g

if [ ! -n "$ncpus" ]; then
    ncpus=40
    echo [`date`]INFO: $ncpus cpus used!
fi


#./gatk+beagle5.sh allbam.list 10

BAM_LIST=$filename

if [ ! -s "${FILE}.raw.vcf.gz" ]; then
    java -Xmx60g -Djava.io.tmpdir=${TEMP_DIR} -jar $GATK \
        -T UnifiedGenotyper \
        -R ${REF_FASTA_FILE} \
        -I $BAM_LIST \
        -nt 24 \
        -nct 2 \
        -L $chr \
        --sample_ploidy 2 \
        -allowPotentiallyMisencodedQuals \
        -out_mode EMIT_VARIANTS_ONLY \
        --genotype_likelihoods_model SNP \
        -o ${FILE}.chr$chr.raw.vcf.gz \
        &> ${FILE}.chr$chr.raw.vcf.gz.log
fi

if [ ! -s ${FILE}.chr$chr.filter.vcf.gz ]; then
    echo [`date`]INFO: "=====" filterSNP
    java -Xmx8g -jar $GATK \
		-T VariantFiltration \
		-R ${REF_FASTA_FILE} \
		-V ${FILE}.chr$chr.raw.vcf.gz \
		-o ${FILE}.chr$chr.filter.vcf.gz \
		--filterExpression "QD < 2.0" --filterName "QDFilter" \
		--filterExpression "MQ < 40.0" --filterName "MQFilter" \
		--filterExpression "QUAL < 50.0 && QUAL >= 30.0" --filterName "LowQual" \
		--filterExpression "QUAL < 30.0" --filterName "VeryLowQual" \
		--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" \
		--filterExpression "FS > 60.0" --filterName "FSFilter" \
		--filterExpression "MQRankSum < -12.5" --filterName "MQRFilter" \
		--filterExpression "ReadPosRankSum < -8.0" --filterName "RPRSFilter" \
		--filterExpression "SOR > 4.0" --filterName "SORFilter" \
		--filterExpression "HaplotypeScore > 13.0" --filterName "HSFilter" \
		-cluster 3 \
		-window 10 \
		--missingValuesInExpressionsShouldEvaluateAsFailing \
		--logging_level ERROR
else
    echo [`date`]INFO: ">>>>>" [${FILE}.filtered.vcf.gz] exists! continue to next step...
fi


if [ ! -s "${FILE}.chr$chr.pass.pos.txt" ]; then
    bcftools view -m2 -M2 -v snps  ${FILE}.chr$chr.filter.vcf.gz  | grep -v "Filter" | grep -v "LowQual" | grep -v "HARD_TO_VALIDATE"  > ${FILE}.chr$chr.pass.vcf
fi

if [ ! -s "${FILE}.chr$chr.pass.vcf.gz" ]; then	
	timer_start=`date "+%Y-%m-%d %H:%M:%S"`	
	java -Xmx100g -XX:-UseGCOverheadLimit -jar ~/software/beagle.22Jul22.46e.jar gt=${FILE}.chr$chr.pass.vcf out=${FILE}.chr$chr.phased nthreads=24
	
	timer_end=`date "+%Y-%m-%d %H:%M:%S"`
	duration=`echo $(($(date +%s -d "${timer_end}") - $(date +%s -d "${timer_start}"))) | awk '{t=split("60 s 60 m 24 h 999 d",a);for(n=1;n<t;n+=2){if($1==0)break;s=$1%a[n]a[n+1]s;$1=int($1/a[n])}print s}'`
	echo "time： $duration" >time.txt		
	gzip ${FILE}.chr$chr.pass.vcf
	
fi	

python2 ./imputeAccuracyID.py  ${FILE}.chr$chr.phased.vcf.gz ${FILE}.chr$chr.pass.vcf.gz /data/songhl/Seq-500/Datadeal/bam/AllCallSNP/${FILE}.chr$chr.pass.vcf.gz /data/songhl/Seq-500/Datadeal/bam/Imputation/ValidationID/id.txt gatk+bg5.txt

exit 0