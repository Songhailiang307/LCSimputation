#!/bin/bash
#PBS -q unlimited
#PBS -l nodes=1:ppn=2
#PBS -N basevar+stitch

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
max_mem=30g

if [ ! -n "$ncpus" ]; then
    ncpus=40
    echo [`date`]INFO: $ncpus cpus used!
fi

#./basevar+stitch.sh allbam.list 10

BAM_LIST=$filename

if [ ! -s "${FILE}.chr$chr.raw.vcf.gz" ]; then

	basevar basetype -R ${REF_FASTA_FILE} \
	--regions $chr \
	--batch-count 50 \
	-L $BAM_LIST \
	--output-vcf ${FILE}.chr$chr.raw.vcf.gz \
	--output-cvg ${FILE}.chr$chr.raw.cvg.tsv.gz \
	--nCPU 40 && echo "** job done **"

rm -r *_temp_*
	
fi



if [ ! -s "${FILE}.chr$chr.filter.vcf.gz" ]; then
	bcftools filter -s LowQual -e '%QUAL<100' \
	-Oz ${FILE}.chr$chr.raw.vcf.gz > ${FILE}.chr$chr.filter.vcf.gz
fi

if [ ! -s "${FILE}.chr$chr.pass.pos.txt" ]; then
    bcftools view -m2 -M2 -v snps  ${FILE}.chr$chr.filter.vcf.gz  | grep -v "Filter" | grep -v "LowQual" | grep -v "HARD_TO_VALIDATE" | bgzip > ${FILE}.chr$chr.pass.vcf.gz
    zcat ${FILE}.chr$chr.pass.vcf.gz| grep "AC=" | awk '{print$1"\t"$2"\t"$4"\t"$5}' > ${FILE}.chr$chr.pass.pos.txt  
fi

#Remove lowercase letters
sed -i '/[a-z]/d' ${FILE}.chr$chr.pass.pos.txt

source ~/.bashrc 
conda activate myrenv
timer_start=`date "+%Y-%m-%d %H:%M:%S"`	
Rscript ./run_STITCH.R $chr $BAM_LIST ${FILE}.chr$chr.pass.pos.txt
timer_end=`date "+%Y-%m-%d %H:%M:%S"`
duration=`echo $(($(date +%s -d "${timer_end}") - $(date +%s -d "${timer_start}"))) | awk '{t=split("60 s 60 m 24 h 999 d",a);for(n=1;n<t;n+=2){if($1==0)break;s=$1%a[n]a[n+1]s;$1=int($1/a[n])}print s}'`
echo "耗时： $duration" >time.txt
conda deactivate

python2 ./imputeAccuracyID.py stitch.$chr.vcf.gz ${FILE}.chr$chr.pass.vcf.gz /data/songhl/Seq-500/Datadeal/bam/AllCallSNP/${FILE}.chr$chr.pass.vcf.gz /data/songhl/Seq-500/Datadeal/bam/Imputation/ValidationID/id.txt base+stitch.txt

	
exit 0