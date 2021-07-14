#BD Rhapsody scRNA-seq mapping manual

#move on to the Rhapsody working directory
cd ~/Rhapsody_analysis

#edit workflow_Rhapsody.sh (specify file name)

#run as follows
nohup sh workflow_Rhapsody.sh in "Rhapsody_analysis" directory.

#parameters of ./shell_scripts/Rhapsody_annotate_noUMI_200cycle_new2.sh (Bowtie2　mapping)

#<Read1> <Read2> : fastq file name. must remove "_L00X_" (lane info)
#Example file name : <samplename>_S<X>_R1_001.fastq.gz <samplename>_S<X>_R2_001.fastq.gz
<samplename> must not contain "_".
<X> is any positive integer value

#t : number of Threads. Must of the half of the theoretical processor number.

#species : Type of organisms. mmu(Mus musculus) or hsa(Homo sapiens).

#adapter : Type of DNA adapter. BDWTA(default, cDNA library) or LibA(old adapter, cDNA library) or hashtag(hashtag), or sampletag(sampletag), Streptavidin(streptavidin)

#BAM : BAM file output flag. 0(default, no BAM output) or 1(enable BAM output)

#Nlane : number of lanes of Novaseq sequencer. If files are L001-L004, specify as "4". If files are L001-L002, specify as "2".

#in_dir : location of input file directory
#out_intermediate : location of intermediate file directory (e.g. concat files, adapter-trimmed files)

#index : index file name for Bowtie2 mapping. GRCm38_101_ENSEMBL(default, for mouse), GRCh38_101_ENSEMBL(human), humanHashtag14(hashtag-human), humanSampleTag12(sampletag-human), mouseHashtag15(hashtag-mouse), mouseSampleTag_BD(sampletag-mouse), Streptavidin_Hashtag(streptavidin).

#data will be outputted into "result" folder.
#do NOT delete "data" folder(temporary store fastq data).



