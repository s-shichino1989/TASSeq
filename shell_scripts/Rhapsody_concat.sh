#!/bin/bash -i

#require cutadapt v2.10
#BD rhapsody WTA proprocessing workflow for ubuntu 16.04, 18.04 or CentOS7, python3 and perl environment
#Written by Shigeyuki Shichino at 20200422


file1=`basename $1 _count_R2.txt | sed -e 's/gathered_//g'`
file2=`basename $2 _count_R2.txt | sed -e 's/gathered_//g'`
file3=`basename $3 _count_R2.txt | sed -e 's/gathered_//g'`
file4=`basename $4 _count_R2.txt | sed -e 's/gathered_//g'`

#pre-count large tidy data

cat $1 $2 $3 $4 | mawk '{sum[$1"\t"$2]+=$3} END {for (name in sum) print name"\t"sum[name]}'> ${file1}_${file2}_${file3}_${file4}_count_R2.txt

sleep 1

rm $1
rm $2
rm $3
rm $4

exit 0
