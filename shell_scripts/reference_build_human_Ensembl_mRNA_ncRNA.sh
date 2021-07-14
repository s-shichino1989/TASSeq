#release-101 build
#build date 20210127
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz


unpigz Homo_sapiens.GRCh38.cdna.all.fa.gz
unpigz Homo_sapiens.GRCh38.ncrna.fa.gz
cat Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.ncrna.fa > Homo_sapiens.GRCh38.rna.all.fa

sleep 1

more Homo_sapiens.GRCh38.rna.all.fa | sed -e 's/^>\([^ \f\n\r\t]*\).*gene_symbol:\([^ \f\n\r\t]*\).*/>\2_\1/g' > Human_RNA.fa

sleep 1

threads=`grep -ce '^processor\s\+:' /proc/cpuinfo`
tmp=`echo $((threads/2))`
bowtie2-build --threads ${tmp} -f Human_RNA.fa GRCh38_101_ENSEMBL

sleep 1

mv GRCh38_101_ENSEMBL.1.bt2 ./index/GRCh38_101_ENSEMBL.1.bt2
mv GRCh38_101_ENSEMBL.2.bt2 ./index/GRCh38_101_ENSEMBL.2.bt2
mv GRCh38_101_ENSEMBL.3.bt2 ./index/GRCh38_101_ENSEMBL.3.bt2
mv GRCh38_101_ENSEMBL.4.bt2 ./index/GRCh38_101_ENSEMBL.4.bt2
mv GRCh38_101_ENSEMBL.rev.1.bt2 ./index/GRCh38_101_ENSEMBL.rev.1.bt2
mv GRCh38_101_ENSEMBL.rev.2.bt2 ./index/GRCh38_101_ENSEMBL.rev.2.bt2

sleep 1

rm Homo_sapiens.GRCh38.cdna.all.fa
rm Homo_sapiens.GRCh38.ncrna.fa
rm Human_RNA.fa
rm Homo_sapiens.GRCh38.rna.all.fa

