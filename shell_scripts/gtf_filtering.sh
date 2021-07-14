#!/bin/sh                                                                       

#require GNU parallel module and seqkit                                         
#wrapper of splitting fastq files                                               

grep -Ff $1 $2 > $1_filtered.gtf
sleep 1

rm $1
