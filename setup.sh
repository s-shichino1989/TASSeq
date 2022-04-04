#!/bin/bash -i
#setup script for analysis environment for BD Rhapsody TAS-Seq data
#required software (for Ubuntu20.04).
#If you use CentOS7 substitute apt to yum and adjust package names)


sudo apt-get update -y
sudo apt-get upgrade -y
sudo apt-get install -y parallel
sudo apt-get install -y wget
sudo apt-get install -y curl pigz
sudo apt-get install -y libncurses5-dev
sudo apt-get install -y zlib1g-dev
sudo apt-get install -y libbz2-dev
sudo apt-get install -y liblzma-dev
sudo apt-get install -y gcc
sudo apt-get install -y libtool texinfo dpkg-dev pkg-config build-essential
sudo apt-get install -y make
sudo apt-get install -y mpich
sudo apt-get install -y python3-pip libssl-dev libpcre3 libpcre3-dev
sudo apt-get install -y nim libopenblas-base
sudo apt-get install -y hdf5-tools hdf5-helpers libhdf5-dev libhdf5-doc libhdf5-serial-dev
sudo apt-get install -y libcurl4-openssl-dev libxml2-dev
sudo apt-get install -y libgmp3-dev
sudo apt-get install -y gdebi
sudo apt-get install -y libgsl-dev default-jre libgeos++-dev libgeos-dev libgeos-doc
sudo apt-get install -y libhts-dev libboost-all-dev libmagick++-dev libboost-dev libboost-all-dev
sudo apt-get install -y python2
sudo apt-get autoremove -y

##set python command as python2 (because hisat2-build depends on python2)
sudo update-alternatives --install /usr/bin/python python /usr/bin/python2.7 1

#change permissions for executable files
sudo chmod 774 ./source_file/split_CellBC_bam
sudo chmod 774 ./source_file/split_CellBC_bam_HISAT2
sudo chmod 774 ./source_file/Lighter/lighter

#clean up remained old R libraries
R_version=`apt list --installed | grep r-base/focal,focal,now | sed -e 's/.*\(4.1.1\).*/\1/g'`

if [ $R_version = 4.1.1 ]; then
echo "R 4.1.1 is already installed"

else
 echo "R 4.1.1 is not installed. remove current R if already installed and re-install clean R 4.1.1."
 sudo apt-get purge -y r-base* r-cran-* r-recommended
 sudo apt-get autoremove -y
 sudo rm -R /usr/local/lib/R/site-library
 #install R 4.1.1
 sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
 sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
 sudo apt-get update -y
 sudo apt-get upgrade -y
 sudo apt-get autoremove -y
 sudo apt-get install -y r-base r-base-dev r-recommended
fi

#installed from source and add PATH to local environment

#samtools
if type "samtools" > /dev/null 2>&1
then
 echo "samtools is already installed"
else
 wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
 tar -jxvf samtools-1.9.tar.bz2
 cd ./samtools-1.9/
 ./configure --prefix=/usr/local/
 make -j 8
 sudo make install
 cd ..
 rm samtools-1.9.tar.bz2
 rm -Rf ./samtools-1.9/
fi

#bowtie2 (add symbolic link)
if type "bowtie2" > /dev/null 2>&1
then
 echo "bowtie2 is already installed"
else
 wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-linux-x86_64.zip
 unzip bowtie2-2.4.2-linux-x86_64.zip
 sudo mv bowtie2-2.4.2-linux-x86_64 /usr/local/bowtie2-2.4.2
 sudo ln -fs /usr/local/bowtie2-2.4.2/bowtie2 /usr/bin
 sudo ln -fs /usr/local/bowtie2-2.4.2/bowtie2-build /usr/bin
 rm bowtie2-2.4.2-linux-x86_64.zip
fi

#FastQC (add symbolic link)
if type "fastqc" > /dev/null 2>&1
then
 echo "fastqc is already installed"
else
 wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
 unzip fastqc_v0.11.9.zip
 sudo chmod 774 ./FastQC/fastqc
 sudo mv FastQC /usr/local/FastQC
 sudo ln -fs /usr/local/FastQC/fastqc /usr/bin
 rm fastqc_v0.11.9.zip
fi

#seqkit (add symbolic link)
if type "seqkit" > /dev/null 2>&1
then
 echo "seqkit is already installed"
else
 wget https://github.com/shenwei356/seqkit/releases/download/v0.15.0/seqkit_linux_amd64.tar.gz
 tar zxvf seqkit_linux_amd64.tar.gz
 sudo mv seqkit /usr/local/seqkit
 sudo ln -fs /usr/local/seqkit /usr/bin
 rm seqkit_linux_amd64.tar.gz
fi

#fastp (add symbolic link)
if type "fastp" > /dev/null 2>&1
then
 echo "fastp is already installed"
else
 wget http://opengene.org/fastp/fastp
 sudo chmod a+x ./fastp
 sudo mv fastp /usr/local/fastp
 sudo ln -fs /usr/local/fastp /usr/bin
 rm fastp
fi

#fftw3
wget http://www.fftw.org/fftw-3.3.9.tar.gz
tar zxvf fftw-3.3.9.tar.gz
cd fftw-3.3.9
./configure --prefix=/usr/local/fftw3 CC=gcc MPICC=mpicc F77=gfortran --enable-mpi --enable-threads --enable-shared --enable-static
make 
sudo make install
cd ..
rm fftw-3.3.9.tar.gz
rm -R fftw-3.3.9
sudo apt-get install -y libfftw3-dev

#hisat2 (add symbolic link)
if type "hisat2" > /dev/null 2>&1
then
 echo "hisat2 is already installed"
else
 curl -s https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download > hisat2-2.2.1-Linux_x86_64.zip
 unzip hisat2-2.2.1-Linux_x86_64.zip
 sudo mv hisat2-2.2.1 /usr/local/hisat2-2.2.1
 sudo ln -fs /usr/local/hisat2-2.2.1/hisat2 /usr/bin
 sudo ln -fs /usr/local/hisat2-2.2.1/hisat2-build /usr/bin
 sudo ln -fs /usr/local/hisat2-2.2.1/hisat2_extract_exons.py /usr/bin
 sudo ln -fs /usr/local/hisat2-2.2.1/hisat2_extract_snps_haplotypes_UCSC.py /usr/bin
 sudo ln -fs /usr/local/hisat2-2.2.1/hisat2_extract_snps_haplotypes_VCF.py /usr/bin
 sudo ln -fs /usr/local/hisat2-2.2.1/hisat2_extract_splice_sites.py /usr/bin
 rm hisat2-2.2.1-Linux_x86_64.zip
fi

#subread/featurecounts
if type "featureCounts" > /dev/null 2>&1
then
echo "featureCounts is already installed"
else
 wget https://sourceforge.net/projects/subread/files/subread-2.0.2/subread-2.0.2-Linux-x86_64.tar.gz
 tar xf subread-2.0.2-Linux-x86_64.tar.gz
 sudo mv subread-2.0.2-Linux-x86_64 /usr/local/subread-2.0.2-Linux-x86_64
 sudo ln -fs /usr/local/subread-2.0.2-Linux-x86_64/bin/featureCounts /usr/bin
 sudo ln -fs /usr/local/subread-2.0.2-Linux-x86_64/bin/exactSNP /usr/bin
 sudo ln -fs /usr/local/subread-2.0.2-Linux-x86_64/bin/subindel /usr/bin
 sudo ln -fs /usr/local/subread-2.0.2-Linux-x86_64/bin/subjunc /usr/bin
 sudo ln -fs /usr/local/subread-2.0.2-Linux-x86_64/bin/sublong /usr/bin
 rm subread-2.0.2-Linux-x86_64.tar.gz
fi

#install rstudio
if type "rstudio" > /dev/null 2>&1
then
 echo "rstudio is already installed"
else
 wget https://download1.rstudio.org/desktop/bionic/amd64/rstudio-2021.09.0-351-amd64.deb
 sudo gdebi -n rstudio-2021.09.0-351-amd64.deb
 rm rstudio-2021.09.0-351-amd64.deb
fi

#required python3 modules (pip3)
#python3.7 or higher
python3 -m pip install --user --upgrade pip
python3 -m pip install --user Cython wheel cmake
python3 -m pip install --user argparse regex Cython wheel cmake pysam==0.15.4 argparse Levenshtein numpy==1.21.0 matplotlib==3.5.1 umap-learn matplotlib pandas subprocess32 htseq==0.12.4 fitsne velocyto scVelo numba==0.52.0

python3 -m pip install cutadapt

python3 -m pip install ./Rhapsody_python/

#add pip path
echo 'export PATH="$PATH:/home/$(whoami)/.local/bin"' >> ~/.bashrc

# install pandoc
if type "pandoc" > /dev/null 2>&1
then
 echo "pandoc is already installed"
else
 wget https://github.com/jgm/pandoc/releases/download/2.11.4/pandoc-2.11.4-1-amd64.deb
 sudo gdebi -n pandoc-2.11.4-1-amd64.deb
 rm pandoc-2.11.4-1-amd64.deb
fi

#download SDMTools
wget https://cran.r-project.org/src/contrib/Archive/SDMTools/SDMTools_1.1-221.2.tar.gz
mv SDMTools_1.1-221.2.tar.gz ./source_file/SDMTools_1.1-221.2.tar.gz


#install required R packages
Rscript ./Rscripts/setup.R

#find out python3 path and modify source R scripts for appropriate PATH
tmp=`which python3`
cat ./Rscripts/library_source_Seurat_template.R | sed -e 's/path_to_python3/${tmp}/g' > ./Rscripts/library_source_Seurat.R

#create bowtie2 index for human and mouse

sh ./shell_scripts/reference_build_mouse_Ensembl_mRNA_ncRNA.sh
sh ./shell_scripts/reference_build_human_Ensembl_mRNA_ncRNA.sh
sh ./shell_scripts/reference_build_cellHashing.sh
sh ./shell_scripts/reference_build_rat_Ensembl_mRNA_ncRNA.sh

#create filtered gtf files
sh ./shell_scripts/Ensembl101_gtf_human_filtering.sh
sh ./shell_scripts/Ensembl101_gtf_mouse_filtering.sh
sh ./shell_scripts/Ensembl104_gtf_rat_filtering.sh

#create hisat2 index files
sh ./shell_scripts/reference_build_Hisat2_mouse_Ensembl.sh
sh ./shell_scripts/reference_build_Hisat2_human_Ensembl.sh
sh ./shell_scripts/reference_build_Hisat2_rat_Ensembl.sh

# install mixcr
if type "mixcr" > /dev/null 2>&1
then
 echo "mixcr is already installed"
else
 wget https://github.com/milaboratory/mixcr/releases/download/v3.0.13/mixcr-3.0.13.zip
 unzip mixcr-3.0.13.zip
 sudo mv mixcr-3.0.13 /usr/local/mixcr-3.0.13
 sudo ln -fs /usr/local/mixcr-3.0.13/mixcr /usr/bin
 rm mixcr-3.0.13.zip
fi

#download imgt libraries for mixcr alignment
wget https://github.com/repseqio/library-imgt/releases/download/v6/imgt.202037-1.sv6.json.gz
unpigz imgt.202037-1.sv6.json.gz
sudo mv imgt.202037-1.sv6.json /usr/local/mixcr-3.0.13/libraries/imgt.202037-1.sv6.json

