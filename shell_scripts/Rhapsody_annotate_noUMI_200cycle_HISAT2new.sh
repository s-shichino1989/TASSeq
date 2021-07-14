#!/bin/bash -i

#require GNU parallel module, pigz, cutadapt, seqkit, HISAT2 and python3.7
#BD Rhapsody WTA proprocessing workflow for ubuntu 20.04, python3.7 environment
#Written by Shigeyuki Shichino at 20210704
#HISAT2 mapping & velocyto preeprocessing workflow

CMDNAME=`basename $0`
if [ $# -ne 10 ]; then
 echo "Usage: sh $CMDNAME sh ./shell_scripts/Rhapsody_annotate_noUMI_200cycle_HISAT2new.sh <samplename>_S<X>_R1_001.fastq.gz <samplename>_S<X>_R2_001.fastq.gz t=<threads: number of CPU threads.> species=<species: hsa or mmu> adapter=<adapter: BDWTA or LibA> BAM=<bam: 0 (not created) or 1(created)> Nlane=<Nlane: number of Illumina lanes: 2 or 4> in_dir=<input data-stored directory> index=<HISAT2 index path> gtf=<gzipped gtf file path>"
 echo "fastq.gz files are stored in <input data-stored directory>, read1 is cell barcode reads and read2 is mRNA reads. Arguments must be fullfilled and ordered correctly." 1>&10
 exit 1
fi

file1=`basename $1 .fastq.gz`
file2=`basename $2 .fastq.gz`
file3=`basename $1 .fastq.gz | sed -e 's/_R1/_L001_R1/g'`
file4=`basename $2 .fastq.gz | sed -e 's/_R2/_L001_R2/g'`
file5=`basename $1 .fastq.gz | sed -e 's/_R1/_L002_R1/g'`
file6=`basename $2 .fastq.gz | sed -e 's/_R2/_L002_R2/g'`
file7=`basename $1 .fastq.gz | sed -e 's/_R1/_L003_R1/g'`
file8=`basename $2 .fastq.gz | sed -e 's/_R2/_L003_R2/g'`
file9=`basename $1 .fastq.gz | sed -e 's/_R1/_L004_R1/g'`
file10=`basename $2 .fastq.gz | sed -e 's/_R2/_L004_R2/g'`
fastqcArg1=`basename $2 .fastq.gz | sed -e 's/_R2_001/_L001_R2_/g'`
fastqcArg2=`basename $2 .fastq.gz | sed -e 's/_R2_001/_L001_R1_/g'`
samplename=`basename $2 | sed -e 's/_\S.*_R2_001.fastq.gz//g'`
sampledate=`basename $2 .fastq.gz | sed -e 's/_\S.*_R2_001//g'`
threads=`echo $3 | sed -e 's/t\=//g'`
species=`echo $4 | sed -e 's/species\=//g'`
adapter=`echo $5 | sed -e 's/adapter\=//g'`
bam=`echo $6 | sed -e 's/BAM\=//g'`
Nlane=`echo $7 | sed -e 's/Nlane\=//g'`
in_dir=`echo $8 | sed -e 's/in_dir\=//g'`
index=`echo ${9} | sed -e 's/index\=//g'`
gtf=`echo ${10} | sed -e 's/gtf\=//g'`
gtf_extract=`echo ${10} | sed -e 's/gtf\=//g' | sed -e 's/.gz//g'`

index1=`basename $index`

if [ $Nlane -eq 2 -o $Nlane -eq 4 ]; then
continue
else
echo "Acceptable Nlanes are 2 or 4, stop."
exit 1
fi

if [ $species = hsa -o $species = mmu -o $species = macaca -o $species = rat ]; then
continue
else
echo "Acceptable species are hsa or mmu or EGFP or tdTomato, stop."
exit 1
fi

if [ $((${threads} % 2)) = 1 ]; then
echo "threads number must be odd number, stop."
exit 1
fi

if [ $bam -eq 0 -o $bam -eq 1 ]; then
continue
else
echo "Acceptable BAM flags are 0(no BAM output) or 1(create BAM output), stop."
exit 1
fi

#load required modules
#module load java
#module load parallel
#module load seqkit
#module load Python3


#cutadapt
start_time=`date +%s`
 echo "Start preprocessing BD rhapsody WTA reads"
 echo "Copy and concat fastq files..."
 echo `date '+%y/%m/%d %H:%M:%S'`

if [ $Nlane = 2 ]; then 
cp ${in_dir}${file3}.fastq.gz ./data/$file3.fastq.gz
cp ${in_dir}${file4}.fastq.gz ./data/$file4.fastq.gz
cp ${in_dir}${file5}.fastq.gz ./data/$file5.fastq.gz
cp ${in_dir}${file6}.fastq.gz ./data/$file6.fastq.gz
cat ./data/$file3.fastq.gz ./data/$file5.fastq.gz > ./data/$file1.fastq.gz
cp ./data/$file1.fastq.gz ${out_intermediate}${file1}.fastq.gz
sleep 1
rm ./data/$file1.fastq.gz
cat ./data/$file4.fastq.gz ./data/$file6.fastq.gz > ./data/$file2.fastq.gz
cp ./data/$file2.fastq.gz ${out_intermediate}${file2}.fastq.gz
sleep 1
rm ./data/$file2.fastq.gz

elif [ $Nlane = 4 ]; then
cp ${in_dir}${file3}.fastq.gz ./data/$file3.fastq.gz
cp ${in_dir}${file4}.fastq.gz ./data/$file4.fastq.gz
cp ${in_dir}${file5}.fastq.gz ./data/$file5.fastq.gz
cp ${in_dir}${file6}.fastq.gz ./data/$file6.fastq.gz
cp ${in_dir}${file7}.fastq.gz ./data/$file7.fastq.gz
cp ${in_dir}${file8}.fastq.gz ./data/$file8.fastq.gz
cp ${in_dir}${file9}.fastq.gz ./data/$file9.fastq.gz
cp ${in_dir}${file10}.fastq.gz ./data/$file10.fastq.gz
#cat ./data/$file3.fastq.gz ./data/$file5.fastq.gz ./data/$file7.fastq.gz ./data/$file9.fastq.gz > ./data/$file1.fastq.gz
#cp ./data/$file1.fastq.gz ${out_intermediate}${file1}.fastq.gz
sleep 1
#rm ./data/$file1.fastq.gz
#cat ./data/$file4.fastq.gz ./data/$file6.fastq.gz ./data/$file8.fastq.gz ./data/$file10.fastq.gz > ./data/$file2.fastq.gz
#cp ./data/$file2.fastq.gz ${out_intermediate}${file2}.fastq.gz
sleep 1
#rm ./data/$file2.fastq.gz

else
 echo "require valid lane number, acceptable is 2 or 4"
 exit 1
fi

echo "Remove 5prime Tagging Adapters from R2 reads..."

if [ $Nlane = 2 ]; then 
ls ./data/${samplename}*_001.fastq.gz | cat | sort | parallel -P 2 -N2 -a \
- 'sh ./shell_scripts/Rhapsody_cutadapt_new.sh' {1} {2} $threads $adapter
sleep 2

#cat ${file3}_trim.fastq.gz ${file5}_trim.fastq.gz > ./data/${file1}_trim.fastq.gz
sleep 1
#cp ./data/${file1}_trim.fastq.gz ${out_intermediate}${file1}_trim.fastq.gz
sleep 1
#rm ./data/${file1}_trim.fastq.gz

#cat ${file4}_trim.fastq.gz ${file6}_trim.fastq.gz > ./data/${file2}_trim.fastq.gz
sleep 1
#cp ./data/${file2}_trim.fastq.gz ${out_intermediate}${file2}_trim.fastq.gz
sleep 1
#rm ./data/${file2}_trim.fastq.gz


elif [ $Nlane = 4 ]; then
tmp=`echo $((threads/4))`
ls ./data/${samplename}*_001.fastq.gz | cat | sort | parallel -P 4 -N2 -a \
- 'sh ./shell_scripts/Rhapsody_cutadapt_new.sh' {1} {2} $tmp $adapter
sleep 2

#cat ${file3}_trim.fastq.gz ${file5}_trim.fastq.gz ${file7}_trim.fastq.gz ${file9}_trim.fastq.gz > ./data/${file1}_trim.fastq.gz
sleep 1
#cp ./data/${file1}_trim.fastq.gz ${out_intermediate}${file1}_trim.fastq.gz
sleep 1
#rm ./data/${file1}_trim.fastq.gz

#cat ${file4}_trim.fastq.gz ${file6}_trim.fastq.gz ${file8}_trim.fastq.gz ${file10}_trim.fastq.gz > ./data/${file2}_trim.fastq.gz
sleep 1
#cp ./data/${file2}_trim.fastq.gz ${out_intermediate}${file2}_trim.fastq.gz
sleep 1
#rm ./data/${file2}_trim.fastq.gz

else
 echo "require valid lane number, acceptable is 2 or 4"
 exit 1
fi

Rscript ./Rscripts/cutadaptlog_concatenate.R ${samplename}

end_time=`date +%s`
time_cutadapt=$((end_time - start_time))

#splitfastq
start_time=`date +%s`
#mkdir ./temp/

if [ $Nlane = 2 ]; then

echo "Start Splitting Fastq into ${threads} files by Seqkit..."
echo `date '+%y/%m/%d %H:%M:%S'`
ls ${samplename}*_001_trim.fastq.gz | cat | sort | parallel -P 2 -N2 -a \
- 'sh ./shell_scripts/Rhapsody_seqkit.sh' {1} {2} $threads
sleep 1

rm ${file3}_trim.fastq.gz
rm ${file4}_trim.fastq.gz
rm ${file5}_trim.fastq.gz
rm ${file6}_trim.fastq.gz

elif [ $Nlane = 4 ]; then

tmp=`echo $((threads/2))`
echo "Start Splitting Fastq into ${tmp} files by Seqkit..."
echo `date '+%y/%m/%d %H:%M:%S'`
ls ${samplename}*_001_trim.fastq.gz | cat | sort | parallel -P 4 -N2 -a \
- 'sh ./shell_scripts/Rhapsody_seqkit.sh' {1} {2} $tmp
sleep 1

rm ${file3}_trim.fastq.gz
rm ${file5}_trim.fastq.gz
rm ${file7}_trim.fastq.gz
rm ${file9}_trim.fastq.gz
rm ${file4}_trim.fastq.gz
rm ${file6}_trim.fastq.gz
rm ${file8}_trim.fastq.gz
rm ${file10}_trim.fastq.gz

else
echo "require valid lane number, acceptable is 2 or 4"
exit 1
fi

sleep 2

#change filename
files='*.part*.fastq.gz'

for filepath in $files
do
tempname1=`basename ${filepath} .fastq.gz |sed -e 's/.*_trim.part_//g'`
tempname2=`basename ${filepath} .fastq.gz |sed -e 's/_001_trim.part_.*/_/g'`
mv $filepath ${tempname1}_${tempname2}.fastq.gz

done

sleep 1
end_time=`date +%s`
time_splitfastq=$((end_time - start_time))


#check fastq quality by FastQC by using subset of fastq data
fastqc -t 2 --nogroup -o ./data 001_${fastqcArg1}.fastq.gz 001_${fastqcArg2}.fastq.gz
sleep 1
unzip ./data/001_${fastqcArg1}.fastqc.zip
unzip ./data/001_${fastqcArg2}.fastqc.zip
mv ./data/001_${fastqcArg1}.fastqc/images/per_base_sequence_content.png ./data/${samplename}_per_base_sequence_content_R1.png
mv ./data/001_${fastqcArg2}.fastqc/images/per_base_sequence_content.png ./data/${samplename}_per_base_sequence_content_R2.png
module unload FastQC
rm ./data/001_${fastqcArg1}.fastqc.zip
rm ./data/001_${fastqcArg2}.fastqc.zip
rm ./data/001_${fastqcArg1}.fastqc.html
rm ./data/001_${fastqcArg2}.fastqc.html
rm -Rf ./data/001_${fastqcArg1}.fastqc
rm -Rf ./data/001_${fastqcArg2}.fastqc


#annotateR1 in parallel
start_time=`date +%s`
echo "Start R1 Annotation..."
echo `date '+%y/%m/%d %H:%M:%S'`

ls *${samplename}*_R1_.fastq.gz | cat | sort | parallel -P ${threads} -a \
- 'python3 ./Rhapsody_python/Rhapsody_python/AnnotateR1_gzip.py' --R1 {}
sleep 1

echo "End R1 Annotation..."
end_time=`date +%s`
time_annotateR1=$((end_time - start_time))


#AnnotateR2
start_time=`date +%s`
echo "Start R2 mapping"
echo `date '+%y/%m/%d %H:%M:%S'`

tmp=`echo $((threads/4))`
ls *${samplename}*_R2_.fastq.gz | cat | sort | parallel -P ${tmp} -a \
- 'sh ./shell_scripts/Rhapsody_hisat2.sh' {} $species $index
sleep 1

Rscript ./Rscripts/mappinglog_concatenate.R ${samplename}
sleep 1
rm *_samtools.log

echo "End R2 mapping by HISAT2"

echo "Start cell BC and UMI assignment to R2 Reads and R2 annotation by featurecounts..."
echo `date '+%y/%m/%d %H:%M:%S'`

unpigz $gtf

sleep 1

echo "Start adding R1 information and R2 annotation and by featureCounts"


ls *${samplename}*_mapping_R2.BAM | cat | sort | parallel -P ${threads} -a \
- 'sh ./shell_scripts/Rhapsody_featureCounts_addtoBam.sh' {} $species $bam $gtf_extract
sleep 1

pigz -p 8 $gtf_extract


Rscript ./Rscripts/mappinglog_concatenate.R ${samplename}
sleep 1
rm *_samtools.log

echo "End R2 annotation"


echo "concatenate bam files..."
ls *${samplename}*_Annotated_mapping_R2.BAM | cat | xargs samtools merge -1 -@${threads} ${samplename}.BAM
sleep 1

rm *${samplename}*_Annotated_mapping_R2.BAM

samtools sort -@24 -m 8G ${samplename}.BAM > ${samplename}.bam
#cp ${samplename}.bam ./result/${samplename}.bam

rm -f ${samplename}.BAM

echo "concatenate annotated files..."
echo `date '+%y/%m/%d %H:%M:%S'`
cat *${samplename}*_count_R2.txt | mawk '{sum[$1"\t"$2]+=$3} END {for (name in sum) print name"\t"sum[name]}'> ./result/${samplename}_final_count.txt

sleep 1
pigz -f ./result/${samplename}_final_count.txt
rm -f *_R2.txt
end_time=`date +%s`
time_annotateR2_and_AddtoSam=$((end_time - start_time))

echo "generate gene-cells martix tables..."
start_time=`date +%s`
sleep 1

#summarize count table (no UMI compression)
Rscript -e "rmarkdown::render('./Rscripts/BDWTA_matrixCreate_HISAT2.Rmd', output_file = '${samplename}_mapping_report.html', output_dir='./result/${samplename}_results/', params=list(samplename='${samplename}', species='${species}', reference_name='${index1}'))"
rm -R ./result/${samplename}_results/${samplename}_mapping_report_files

sleep 1
rm ./result/*${samplename}*.log


#create cell barcode-splitted BAM files and loom file for velocyto analysis
if [ $bam = 0 ]; then
rm ${samplename}.csv

else 
Rscript ./Rscripts/CSVsplit.R ${samplename}.csv ${threads}
rm ${samplename}.csv

ls ${samplename}*.csv | parallel -P ${threads} -a \
- './source_file/split_CellBC_bam_HISAT2' {} ${samplename}.bam ${samplename}

mkdir ./result/${samplename}_BAM

rm -f ${samplename}.bam

sleep 1

ls *.bam | parallel -P ${threads} -a \
- 'sh ./shell_scripts/filemove.sh' {} ./result/${samplename}_BAM

sleep 1

rm -f *${samplename}*.bam
rm -f ${samplename}*.csv 

#make loom file for velocyto analysis
sh ./shell_scripts/velocyto_script.sh ${gtf} ./result/${samplename}_BAM ${threads}

tar -cf ./result/${samplename}_BAM.tar ./result/${samplename}_BAM
#pigz -f ./result/${samplename}_BAM.tar

sleep 1

rm -Rf ./result/${samplename}_BAM
rm *.BAM

mv ./result/${samplename}_BAM.tar ./result/${samplename}_results/processed/${samplename}_BAM.tar
md5sum ./result/${samplename}_results/processed/${samplename}_BAM.tar > ./result/${samplename}_results/processed/${samplename}_BAM_md5sum.txt

#move processed files
mv ./result/${samplename}_BAM_velocyto_combined.loom.gz ./result/${samplename}_results/processed/${samplename}_velocyto_combined.loom.gz
md5sum ./result/${samplename}_results/processed/${samplename}_velocyto_combined.loom.gz > ./result/${samplename}_results/processed/${samplename}_velocyto_md5sum.txt
rm -f ./result/${samplename}.bam

fi

#move processed files
mv ./result/${samplename}_final_count.txt.gz ./result/${samplename}_results/processed/${samplename}_final_count.txt.gz

mv ./data/${samplename}_per_base_sequence_content_R1.png ./result/${samplename}_results/plots/${samplename}_per_base_sequence_content_R1.png
mv ./data/${samplename}_per_base_sequence_content_R2.png ./result/${samplename}_results/plots/${samplename}_per_base_sequence_content_R2.png

end_time=`date +%s`
time_reshaping=$((end_time - start_time))

echo "Preprocessing finished." 
echo "Cutadapt ${time_cutadapt} seconds" 
echo "splitfastq ${time_splitfastq} seconds" 
echo "annotateR1 ${time_annotateR1} seconds" 
echo "annotateR2_and_AddtoSam ${time_annotateR2_and_AddtoSam} seconds" 
echo "generating matrix ${time_reshaping} seconds" 
time=$((time_cutadapt + time_splitfastq + time_annotateR1 + time_annotateR2_and_AddtoSam + time_reshaping))
echo "Total running time is ${time} seconds"
echo `date '+%y/%m/%d %H:%M:%S'`
#module purge
exit 0
