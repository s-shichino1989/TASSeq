#!/bin/bash -i

#require GNU parallel module, velocyto.py and loompy

GTF_FILE=`basename $1 .gtf.gz`
BAM_dir=$2
threads=$3

mkdir ${BAM_dir}_results
sleep 1

cp $1 ${BAM_dir}_${GTF_FILE}.gtf.gz

unpigz ${BAM_dir}_${GTF_FILE}.gtf.gz

#module load parallel
#module load Python3

ls ${BAM_dir}/*.bam | cat | sort | parallel -P $threads -a \
- 'sh ./shell_scripts/velocyto_internal.sh' {} ${BAM_dir}_results ${BAM_dir}_${GTF_FILE}.gtf $BAM_dir

#concatenate loom files
ls ${BAM_dir}_results/*.loom | cat > ${BAM_dir}_loom_files.txt

sleep 1

python3 ./external_softwares/loompy.combine.py --filelist ${BAM_dir}_loom_files.txt --out ${BAM_dir}_velocyto_combined.loom

sleep 1

#compress loom file
pigz -p 8 ${BAM_dir}_velocyto_combined.loom

rm ${BAM_dir}_loom_files.txt
rm ${BAM_dir}_${GTF_FILE}.gtf
#rm -Rf ${BAM_dir}_results




