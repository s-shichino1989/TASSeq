#!/bin/bash -i

#BD rhapsody WTA proprocessing workflow, Hisat2 mapping module
#Written by Shigeyuki Shichino at 20200703

CMDNAME=`basename $0`
if [ $# -ne 3 ]; then
 echo "Usage: sh $CMDNAME Read2.fastq.gz <species>"
 echo "read2 is mRNA reads"
 echo "example : sh $CMDNAME day00-1_S1_R2_001.fastq.gz hsa index" 1>&3
 exit 1
fi

file1=`basename $1 _R2_.fastq.gz`
species=`echo $2`
index=`echo $3`

#load required modules
#module load hisat2
#module load samtools

#mapping reads with hisat2, remove secondary alignments

if [ $species = mmu ]; then
hisat2 -q -p 6 --rna-strandness F --very-sensitive --seed 656565 --reorder --omit-sec-seq --mm \
-x $3 -U $1 | grep -v "0\s\*\s\*\sAS:i:" - | 
samtools view -@2 -Shb - > ${file1}_mapping_R2.BAM
rm $1

elif [ $species = hsa ]; then
hisat2 -q -p 6 --rna-strandness F --very-sensitive --seed 656565 --reorder --omit-sec-seq --mm \
-x $3 -U $1 | grep -v "0\s\*\s\*\sAS:i:" - | 
samtools view -@2 -Shb - > ${file1}_mapping_R2.BAM

rm $1

else
 echo "require valid species, acceptable is mmu or hsa"
 exit 1
fi

samtools flagstat -@4 ${file1}_mapping_R2.BAM | sed '2,4d' - | sed '3,10d' | sed -e 's/ .*//g' > ${file1}_samtools.log

exit 0
