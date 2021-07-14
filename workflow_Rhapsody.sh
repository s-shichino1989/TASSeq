#!/bin/bash -i

#usage : nohup sh workflow_Rhapsody.sh

#If you want to perform annotation/mapping of mRNA reads or hashtag reads or sampletag reads
#sh ./shell_scripts/Rhapsody_annotate_noUMI_200cycle_new2.sh <samplename>_S<X>_R1_001.fastq.gz <samplename>_S<X>_R2_001.fastq.gz t=24 species=mmu adapter=BDWTA BAM=0 Nlane=2 in_dir=/HDD1/data/ out_intermediate=/HDD1/trimmed/ index=./index/GRCm38_101_ENSEMBL

#If you want to perform demultiplexing
#Rscript -e "rmarkdown::render('./Rscripts/demultiplex_only.Rmd', output_file = '<samplename>_demulti_report.html', output_dir='./result/<samplename>_results/', params=list(cDNA_matrix='matrix_inflection_<samplename>.txt.gz', Tag_matrix='Hashtag_top1M_<samplenameTag>.txt.gz', nTags='3', threads='4', DBEC_output='FALSE', species='hsa', min.ave='6', min.diff='5.5', SeqGeq_output='FALSE'))"
