#!/bin/bash -i

#internal function of parallel processing for loom file generation by velocyto analysis

filename=`basename $1 .bam`
result_path=$2
GTF_file=$3
input_path=$4

#module load Python3

velocyto run -@ 1 -c -U -o ${result_path}/${filename} ${input_path}/${filename}.bam ${GTF_file}

sleep 1

mv ${result_path}/${filename}/${filename}*.loom ${result_path}/${filename}.loom 

sleep 1

rm -R ${result_path}/${filename}

