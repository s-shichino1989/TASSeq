#!/bin/bash -i
#setup script for analysis environment for BD Rhapsody TAS-Seq data
#for Ubuntu20.04.
#If you use CentOS, substitute apt to yum and adjust package names.

sudo apt-get update -y
sudo apt-get upgrade -y
sudo apt-get install -y parallel
sudo apt-get install -y wget
sudo apt-get install -y curl
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
sudo apt-get install -y libgsl-dev default-jre
sudo apt-get install -y libhts-dev

#clean up remained old R libraries
R_version=`apt info r-base | grep "Version" | sed -e 's/.*\(3.6.3\).*/\1/g'`

if [ $R_version = 3.6.3 ]; then
echo "R 3.6.3 is already installed"

else
 echo "R 3.6.3 is not installed. remove current R if already installed and re-install clean R 3.6.3."
 sudo apt-get purge -y r-base* r-cran-* r-recommended
 sudo apt-get autoremove -y
 sudo rm -R /usr/local/lib/R/site-library
 #install R 3.6.3 (do NOT install R 4.0)
 sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
 sudo add-apt-repository --remove 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
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
 wget https://github.com/shenwei356/seqkit/releases/download/v0.14.0/seqkit_linux_amd64.tar.gz
 tar zxvf seqkit_linux_amd64.tar.gz
 sudo mv seqkit /usr/local/seqkit
 sudo ln -fs /usr/local/seqkit /usr/bin
 rm seqkit_linux_amd64.tar.gz
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

#install rstudio
if type "rstudio" > /dev/null 2>&1
then
 echo "rstudio is already installed"
else
 wget https://download1.rstudio.org/desktop/bionic/amd64/rstudio-1.4.1103-amd64.deb
 sudo gdebi -n rstudio-1.4.1103-amd64.deb
 rm rstudio-1.4.1103-amd64.deb
fi

#required python3 modules (pip3)
#python3.7 or higher
sudo python3 -m pip install --upgrade pip
python3 -m pip install --user --upgrade pip
sudo python3 -m pip install cutadapt Cython wheel cmake
python3 -m pip install argparse regex cutadapt \
pysam==0.15.4 argparse Levenshtein numpy umap-learn \
matplotlib pandas subprocess32 htseq==0.12.4 fitsne MulticoreTSNE \
velocyto scVelo numba==0.52.0

python3 -m pip install ./Rhapsody_python/

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
sudo Rscript ./Rscripts/setup.R

#find out python3 path and modify source R scripts for appropriate PATH
tmp=`which python3`
cat ./Rscripts/library_source_Seurat_template.R | sed -e 's/path_to_python3/${tmp}/g' > ./Rscripts/library_source_Seurat.R


#create bowtie2 index for human and mouse

sh ./shell_scripts/reference_build_mouse_Ensembl_mRNA_ncRNA.sh
sh ./shell_scripts/reference_build_human_Ensembl_mRNA_ncRNA.sh
sh ./shell_scripts/reference_build_cellHashing.sh


