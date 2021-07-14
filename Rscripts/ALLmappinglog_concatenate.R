#!/usr/bin/env Rscript

fnames = commandArgs(trailingOnly=TRUE)
num_of_data = length(fnames)

res=read.table(fnames[1], header=TRUE, row.names=1, sep="\t", quote="")

for (i in 2:num_of_data){
 tmp = read.table(fnames[i], header=TRUE, row.names=1, sep="\t", quote="")
 res=cbind(res,tmp)
 }

fnames1 = unlist(fnames)
fnames1 = as.character(fnames1)
fnames1 = gsub(".txt", "", fnames1)
fnames1 = gsub("./result/stats_concat_", "", fnames1)
colnames(res)=fnames1
invisible(file.remove(as.character(unlist(fnames))))

write.table(res, "./result/ALLstats.log", row.names=T, sep="\t", quote=F)



