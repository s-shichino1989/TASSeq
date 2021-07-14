#!/bin/bash -i

#BD rhapsody WTA proprocessing workflow, Hisat2 mapping module
#Written by Shigeyuki Shichino at 20200703

CMDNAME=`basename $0`
if [ $# -ne 4 ]; then
 echo "Usage: sh $CMDNAME Read2.fastq.gz <species> <BAM flag> <gtf_file>"
 echo "read2 is mRNA reads" 
 echo "example : sh $CMDNAME day00-1_S1_R2_001.fastq.gz hsa 0 gtf" 1>&4
 exit 1
fi

file1=`basename $1 _mapping_R2.BAM`

species=`echo $2`
gtf=`echo $4`

#module load samtools
#module load Python3

#retain mapping info, cell BC and MI-containing BAM files 

python3 ./Rhapsody_python/Rhapsody_python/AddtoBam_gzip.py --annot-R1 ${file1}_Annotation_R1.csv.gz --bam ${file1}_mapping_R2.BAM

sleep 1

#STAR/HISAT2-mapped reads are annotated by FeatureCounts
featureCounts -T 2 -Q 0 -s 1 -t gene -g gene_name \
-a $gtf --primary -M -O --largestOverlap --fraction -R BAM -o ${file1}_counts.txt ${file1}_Annotated_mapping_R2.BAM

sleep 1
rm ${file1}_counts.txt
rm ${file1}_counts.txt.summary
rm ${file1}_Annotated_mapping_R2.BAM



samtools flagstat -@2 ${file1}_Annotated_mapping_R2.BAM.featureCounts.bam | sed '2,4d' - | sed '3,10d' | sed -e 's/ .*//g' > ${file1}_samtools.log
sleep 1

#remove secondary annotation from bam file
samtools view -@2 -h ${file1}_Annotated_mapping_R2.BAM.featureCounts.bam | sed -e 's/\(XT:Z:[^,]*\),[^\t\s\n]*/\1/g' - | samtools view -@4 -Sb | samtools sort -@2 -m 2G - > ${file1}_Annotated_mapping_R2.BAM

rm ${file1}_Annotated_mapping_R2.BAM.featureCounts.bam


#count by gene features and cell BCs

samtools view -@4 ${file1}_Annotated_mapping_R2.BAM | grep -v -e '^@\|^#\|^$' - | grep -v 'CB:Z:0' | grep 'XS:Z:Assigned' | mawk '{print $NF"\t"$(NF-4)}' | grep -v '^\*' | sed -e 's/XT:Z://g' -e 's/,.*\t/\t/g' -e 's/CB:Z://g' | mawk '{count[$0]++} END {for (name in count) print name"\t"count[name]}' > gathered_${file1}_count_R2.txt
sleep 1


exit 0
