#!/usr/bin/env Rscript
#setup R environment

install.packages(c("data.table", "dplyr", "ggplot2", 
"stringr", "tibble", "tidyr", "mcgv", "voxel", "ggplotify","foreach", "doParallel", "withr",
"Matrix", "mclust", "recommenderlab", "MASS", "future.apply", "future", "doFuture", "reticulate", "viridis", "devtools", "factoextra", "flashClust", "cowplot", "XML", "RCurl"), lib="/usr/local/lib/R/site-library")

#install bioconductor packages for R 3.6.3
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version="3.10")
BiocManager::install(c("SingleCellExperiment","SummarizedExperiment", "GenomicRanges", "TCC", "impute", "preprocessCore", "GO.db", "AnnotationDbi", "TCC", "BiocParallel", "AnnotationDbi", "biomaRt"), lib="/usr/local/lib/R/site-library")

BiocManager::install(c("DropletUtils","flowTrans", "slingshot", "clusterExperiment",
"monocle", "multtest", "clusterProfiler"), lib="/usr/local/lib/R/site-library")

library(devtools)
library(withr)

withr::with_libpaths(new="/usr/local/lib/R/site-library", install_version("RcppAnnoy", version = "0.0.16")) #need downgrade RcppAnnoy for install BiocNeighbors
BiocManager::install('BiocNeighbors', type = 'source', lib="/usr/local/lib/R/site-library")
BiocManager::install(c("SingleR"), version="3.10", lib="/usr/local/lib/R/site-library")
#if you use CentOS7 you must install hdf5 manually
#wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.20/src/hdf5-1.8.20.tar.gz
#unpigz hdf5-1.8.20.tar.gz
#tar xf hdf5-1.8.20.tar.gz
#cd hdf5-1.8.20
#./configure --enable-fortran --prefix=/usr/local/hdf5 --enable-cxx --with-pthread=/usr/include/,/usr/lib/x86_64-linux-gnu
#make -j 24
#make check
#sudo make install
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/hdf5/lib
#install.packages("hdf5r", configure.args="--with-hdf5=/usr/local/hdf5/bin/h5cc")
withr::with_libpaths(new="/usr/local/lib/R/site-library", install_url("https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_1.3.0.tar.gz"))

install.packages("WGCNA", "parallelDist", "ggsci", lib="/usr/local/lib/R/site-library")

install.packages("./source_file/SDMTools_1.1-221.2.tar.gz", type="source", repos=NULL, lib="/usr/local/lib/R/site-library")

#install Seurat v2.3.4
source("https://z.umn.edu/archived-seurat")

#install modified tradeSeq package
BiocManager::install("tradeSeq")
install.packages("./source_file/rDBEC_0.2.3.tar.gz", type="source", repos=NULL, lib="/usr/local/lib/R/site-library")

withr::with_libpaths(new="/usr/local/lib/R/site-library", install_github('barkasn/fastSave', lib="/usr/local/lib/R/site-library"))
install.packages(c("rmdformats", "DT", "ggpubr"), lib="/usr/local/lib/R/site-library")
withr::with_libpaths(new="/usr/local/lib/R/site-library", install_github("LTLA/celldex"))

install_github("velocyto-team/velocyto.R")



