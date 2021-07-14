#!/usr/bin/env Rscript
#extract mapping statistics from SAM files

fnames = dir(pattern = "*samtools.log")
num_of_data = length(fnames)
filename = commandArgs(trailingOnly=TRUE)[1] 

res=read.table(fnames[1], header=FALSE, sep="\t", quote="")

for (i in 2:num_of_data){
 tmp = read.table(fnames[i], header=FALSE, sep="\t", quote="")
 res=cbind(res,tmp)
 }

tmp = rowSums(res)
res=cbind(res,tmp)
tmp = res[2,]/res[1,]*100
res = rbind(res,tmp)

fnames = gsub("\\..+$", "", fnames) # remove ".txt" from file names
fnames = gsub("_samtools", "", fnames) # remove "_samtools" from file names
rownames(res)=c("total_reads", "mapped_reads", "mapping_rate")

file.name = sprintf("./result/%s_mapping.log", filename)
write.table(res[,c(num_of_data, num_of_data+1)], file.name, row.names=T, sep="\t", quote=F)



