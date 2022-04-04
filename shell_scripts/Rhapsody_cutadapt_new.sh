#!/bin/bash -i

#require cutadapt v2.10
#BD rhapsody WTA proprocessing workflow for ubuntu 16.04, 18.04 or CentOS7, python3 and perl environment
#Written by Shigeyuki Shichino at 20200422

CMDNAME=`basename $0`
if [ $# -ne 4 ]; then
 echo "Usage: sh $CMDNAME Read1.fastq.gz Read2.fastq.gz t<num_threads> <adapter>"
 echo "fastq.gz files are stored in ./data folder, read1 is cell barcode reads and read2 is mRNA reads" 
 echo "example : sh $CMDNAME day00-1_S1_R1_001.fastq.gz day00-1_S1_R2_001.fastq.gz 16 BDWTA" 1>&4
 exit 1
fi

file1=`basename $1 .fastq.gz`
file2=`basename $2 .fastq.gz`
threads=`echo $3 `
adapter=`echo $4 `
samplename=`basename $2 | sed -e 's/_R2_001.fastq.gz//g'`

#load required modules

#trimming adapters by cutadapt v2.10
if [ $adapter = "LibA" ]; then

cutadapt --cores=$threads -u -1 -U -1 -q 20 --trim-n --report=minimal -Z -n 3 \
-G "XCTATGCGCCTTGCCAGCCCGCTCAGACT;max_error_rate=0.1" -G "G{100};min_overlap=8;max_error_rate=0.2" \
-A "A{10};min_overlap=6;max_error_rate=0.2" -m 60:30 \
-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log &
pid=$!
wait $pid

rm ./data/${file1}.fastq.gz
rm ./data/${file2}.fastq.gz

elif [ $adapter = "BDWTA" ]; then

cutadapt --cores=$threads -u -1 -U -1 -q 20 --trim-n --report=minimal -Z -n 3 \
-G "XAAGCAGTGGTATCAACGCAGA;max_error_rate=0.1" -G "G{100};min_overlap=8;max_error_rate=0.2" \
-A "A{10};min_overlap=6;max_error_rate=0.2" -m 60:30 \
-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log

#cutadapt --cores=$threads -u -84 -U -1 -q 20 --trim-n --report=minimal -Z -n 3 \
#-G "XAAGCAGTGGTATCAACGCAGA;max_error_rate=0.1" -G "G{100};min_overlap=8;max_error_rate=0.2" \
#-A "A{10};min_overlap=6;max_error_rate=0.2" -m 60:30 \
#-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
#./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log



rm ./data/${file1}.fastq.gz
rm ./data/${file2}.fastq.gz

elif [ $adapter = "targeted" ]; then

cutadapt --cores=$threads -u -1 -U -1 -q 20 --report=minimal -Z --trim-n -A "A{10};min_overlap=6;max_error_rate=0.2" -m 60:30 \
-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log

rm ./data/${file1}.fastq.gz
rm ./data/${file2}.fastq.gz


elif [ $adapter = "hashtag" ]; then

cutadapt --cores=$threads -U -100 -q 20 --trim-n -Z --report=minimal -m 60:15 \
-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log

sleep 1
rm ./data/${file1}.fastq.gz
rm ./data/${file2}.fastq.gz

elif [ $adapter = "sampletag" ]; then
cutadapt --cores=$threads -u -1 -U -90 -q 20 --trim-n -Z --report=minimal -G "XGTTGTCAAGATGCTACCGTTCAGAG;min_overlap=20" -m 60:40 \
-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log

#cutadapt --cores=$threads -u -1 -q 20 --trim-n -Z --report=minimal -m 56:40 \
#-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
#./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log

rm ./data/${file1}.fastq.gz
rm ./data/${file2}.fastq.gz


elif [ $adapter = "Streptavidin" ]; then

cutadapt --cores=$threads -U -82 -q 20 --trim-n -Z --report=minimal -G "XTGGCACCCGAGAATTCCA;min_overlap=17" -m 60:15 \
-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log

sleep 1
rm ./data/${file1}.fastq.gz
rm ./data/${file2}.fastq.gz

elif [ $adapter = "Abseq" ]; then

cutadapt --cores=$threads -U 12 -U -58 -q 20 --trim-n -Z --report=minimal -m 60:40 \
-o ${file1}_trim.fastq.gz -p ${file2}_trim.fastq.gz \
./data/${file1}.fastq.gz ./data/${file2}.fastq.gz 1> ${samplename}_cutadapt_summary.log

sleep 1
rm ./data/${file1}.fastq.gz
rm ./data/${file2}.fastq.gz


else
 echo "require valid adapter, acceptable is BDWTA, LibA, targeted, hashtag or sampletag or Streptavidin"
 exit 1
fi

exit 0
