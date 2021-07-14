#!/bin/bash -i

#require GNU parallel module, pigz, cutadapt, bowtie2 and python3
#BD rhapsody WTA proprocessing workflow for ubuntu 16.04, python3 and perl environment
#Written by Shigeyuki Shichino at 20190310

CMDNAME=`basename $0`
if [ $# -ne 10 ]; then
 echo "fastq.gz files are stored in <input data-stored directory>, read1 is cell barcode reads and read2 is mRNA reads. Arguments must be fullfilled and ordered correctly."
 echo "Usage: sh $CMDNAME sh ./shell_scripts/Rhapsody_annotate_noUMI_200cycle_new2.sh <samplenmae2>_S<X>_R1_001.fastq.gz <samplename2>_S<X>_R2_001.fastq.gz t=<threads: number of CPU threads.> species=<species: hsa or mmu> adapter=<adapter: BDWTA or LibA> BAM=<bam: 0 (not created) or 1(created)> Nlane=<Nlane: number of Illumina lanes: 2 or 4> in_dir=<input data-stored directory> out_intermediate=<directory in which concatenated fastq.gz files will be stored> index=<bowtie2 index path>" 1>&10
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
species=`echo $4| sed -e 's/species\=//g'`
adapter=`echo $5| sed -e 's/adapter\=//g'`
bam=`echo $6 | sed -e 's/BAM\=//g'`
Nlane=`echo $7 | sed -e 's/Nlane\=//g'`
in_dir=`echo $8 | sed -e 's/in_dir\=//g'`
out_intermediate=`echo $9 | sed -e 's/out_intermediate\=//g'`
index=`echo ${10} | sed -e 's/index\=//g'`


if [ $Nlane -eq 2 -o $Nlane -eq 4 ]; then
continue
else
echo "Acceptable Nlanes are 2 or 4, stop."
exit 1
fi

if [ $species = hsa -o $species = mmu -o $species = EGFP -o $species = tdTomato -o $species = macaca ]; then
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
module load parallel
module load seqkit
module load Python3

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
cat ./data/$file3.fastq.gz ./data/$file5.fastq.gz ./data/$file7.fastq.gz ./data/$file9.fastq.gz > ./data/$file1.fastq.gz
cp ./data/$file1.fastq.gz ${out_intermediate}${file1}.fastq.gz
sleep 1
rm ./data/$file1.fastq.gz
cat ./data/$file4.fastq.gz ./data/$file6.fastq.gz ./data/$file8.fastq.gz ./data/$file10.fastq.gz > ./data/$file2.fastq.gz
cp ./data/$file2.fastq.gz ${out_intermediate}${file2}.fastq.gz
sleep 1
rm ./data/$file2.fastq.gz

else
 echo "require valid lane number, acceptable is 2 or 4"
 exit 1
fi

echo "Remove 5prime Tagging Adapters from R2 reads..."

if [ $Nlane = 2 ]; then 
tmp=`echo $((threads/2))`
ls ./data/${samplename}*_001.fastq.gz | cat | sort | parallel -P 2 -N2 -a \
- 'sh ./shell_scripts/Rhapsody_cutadapt_new.sh' {1} {2} $threads $adapter
sleep 1

elif [ $Nlane = 4 ]; then
tmp=`echo $((threads/2))`
ls ./data/${samplename}*_001.fastq.gz | cat | sort | parallel -P 4 -N2 -a \
- 'sh ./shell_scripts/Rhapsody_cutadapt_new.sh' {1} {2} $tmp $adapter
sleep 1

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
unzip ./data/001_${fastqcArg1}_fastqc.zip
unzip ./data/001_${fastqcArg2}_fastqc.zip
mv 001_${fastqcArg1}_fastqc/Images/per_base_sequence_content.png ./data/${samplename}_per_base_sequence_content_R1.png
mv 001_${fastqcArg2}_fastqc/Images/per_base_sequence_content.png ./data/${samplename}_per_base_sequence_content_R2.png
module unload FastQC
rm ./data/001_${fastqcArg1}_fastqc.zip
rm ./data/001_${fastqcArg2}_fastqc.zip
rm ./data/001_${fastqcArg1}_fastqc.html
rm ./data/001_${fastqcArg2}_fastqc.html
rm -Rf 001_${fastqcArg1}_fastqc
rm -Rf 001_${fastqcArg2}_fastqc


#annotateR1 in parallel
start_time=`date +%s`
echo "Start R1 Annotation..."
echo `date '+%y/%m/%d %H:%M:%S'`

tmp=`echo $((threads/2))`
ls *${samplename}*_R1_.fastq.gz | cat | sort | parallel -P ${tmp} -a \
- 'python3 ./Rhapsody_python/Rhapsody_python/AnnotateR1_gzip.py' --R1 {}
sleep 1

echo "End R1 Annotation..."
end_time=`date +%s`
time_annotateR1=$((end_time - start_time))

#AnnotateR2 and AddtoSam
start_time=`date +%s`
echo "Start R2 Annotation and add R1 information to mapped R2 files..."
echo `date '+%y/%m/%d %H:%M:%S'`

tmp=`echo $((threads/2))`

ls *${samplename}*_R2_.fastq.gz | cat | sort | parallel -P ${tmp} -a \
- 'sh ./shell_scripts/Rhapsody_bowtie2_addtoBam.sh' {} $species $bam $adapter $index
sleep 1


Rscript ./Rscripts/mappinglog_concatenate.R ${samplename}
sleep 1
rm *_samtools.log

echo "End R1 and R2 Annotation"

if [ $bam = 1 ]; then

echo "concatenate bam files..."
ls *${sampledate}*.BAM | cat | xargs samtools merge -@ ${threads} ${samplename}.BAM 
sleep 1

cp ${samplename}.BAM ${out_intermediate}${samplename}.BAM 

rm *.BAM

else
sleep 1
fi

echo "concatenate annotated files..."
#cat *_count_R2.txt | mawk -F "\t" '{count[$0]++} END {for (name in count) print name"\t"count[name]}' > ./result/${samplename}_final_count.txt 

ls *_count_R2.txt | cat | sort | parallel -P 4 -N4 -a \
- 'sh ./shell_scripts/Rhapsody_concat.sh' {1} {2} {3} {4}
sleep 1

cat *_count_R2.txt | mawk '{sum[$1"\t"$2]+=$3} END {for (name in sum) print name"\t"sum[name]}'> ./result/${samplename}_final_count.txt

sleep 1
pigz -f ./result/${samplename}_final_count.txt
rm *_R2.txt
end_time=`date +%s`
time_annotateR2_and_AddtoSam=$((end_time - start_time))

echo "generate gene-cells martix tables..."
start_time=`date +%s`
sleep 1

#summarize count table (no UMI compression)
#create html report by Rmarkdown package

if [ $adapter != "hashtag" ] && [ $adapter != "sampletag" ] && [ $adapter != "Streptavidin" ]; then

Rscript -e "rmarkdown::render('./Rscripts/BDWTA_matrixCreate.Rmd', output_file = '${samplename}_mapping_report.html', output_dir='./result/${samplename}_results/', params=list(samplename='${samplename}', species='${species}'))"
rm -R ./result/${samplename}_results/${samplename}_report_files

sleep 1
rm ./result/*${samplename}*.log

#move processed files
mv ./result/${samplename}_final_count.txt.gz ./result/${samplename}_results/processed/${samplename}_final_count.txt.gz
mv ./data/${samplename}_per_base_sequence_content_R1.png ./result/${samplename}_results/plots/${samplename}_per_base_sequence_content_R1.png
mv ./data/${samplename}_per_base_sequence_content_R2.png ./result/${samplename}_results/plots/${samplename}_per_base_sequence_content_R2.png

elif [ $adapter = "hashtag" ] || [ $adapter = "sampletag" ] || [ $adapter = "Streptavidin" ]; then

Rscript -e "rmarkdown::render('./Rscripts/BDTag_matrixCreate.Rmd', output_file = '${samplename}_mapping_report.html', output_dir='./result/${samplename}_results/', params=list(samplename='${samplename}', species='${species}', adapter='${adapter}'))"
rm -R ./result/${samplename}_results/${samplename}_report_files

sleep 1
rm ./result/*${samplename}*.log

#move processed files
mv ./result/${samplename}_final_count.txt.gz ./result/${samplename}_results/processed/${samplename}_final_count.txt.gz
mv ./data/${samplename}_per_base_sequence_content_R1.png ./result/${samplename}_results/plots/${samplename}_per_base_sequence_content_R1.png
mv ./data/${samplename}_per_base_sequence_content_R2.png ./result/${samplename}_results/plots/${samplename}_per_base_sequence_content_R2.png

else
sleep 1
fi

end_time=`date +%s`
time_reshaping=$((end_time - start_time))

rm Rplots.pdf


echo "Preprocessing finished." 
echo "Cutadapt ${time_cutadapt} seconds" 
echo "splitfastq ${time_splitfastq} seconds" 
echo "annotateR1 ${time_annotateR1} seconds" 
echo "annotateR2_and_AddtoSam ${time_annotateR2_and_AddtoSam} seconds" 
echo "generating matrix ${time_reshaping} seconds" 
time=$((time_cutadapt + time_splitfastq + time_annotateR1 + time_annotateR2_and_AddtoSam + time_reshaping))
echo "Total running time is ${time} seconds"
echo `date '+%y/%m/%d %H:%M:%S'`
module purge
exit 0
