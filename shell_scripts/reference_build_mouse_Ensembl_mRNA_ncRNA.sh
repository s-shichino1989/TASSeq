#release-101 build
#build date 20210127
wget ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz

unpigz Mus_musculus.GRCm38.cdna.all.fa.gz
unpigz Mus_musculus.GRCm38.ncrna.fa.gz
cat Mus_musculus.GRCm38.cdna.all.fa Mus_musculus.GRCm38.ncrna.fa > Mus_musculus.GRCm38.rna.all.fa

sleep 1

more Mus_musculus.GRCm38.rna.all.fa | sed -e 's/^>\([^ \f\n\r\t]*\).*gene_symbol:\([^ \f\n\r\t]*\).*/>\2_\1/g' > Mouse_RNA.fa

sleep 1

threads=`grep -ce '^processor\s\+:' /proc/cpuinfo`
tmp=`echo $((threads/2))`
bowtie2-build --threads ${tmp} -f Mouse_RNA.fa GRCm38_101_ENSEMBL

sleep 1

mv GRCm38_101_ENSEMBL.1.bt2 ./index/GRCm38_101_ENSEMBL.1.bt2
mv GRCm38_101_ENSEMBL.2.bt2 ./index/GRCm38_101_ENSEMBL.2.bt2
mv GRCm38_101_ENSEMBL.3.bt2 ./index/GRCm38_101_ENSEMBL.3.bt2
mv GRCm38_101_ENSEMBL.4.bt2 ./index/GRCm38_101_ENSEMBL.4.bt2
mv GRCm38_101_ENSEMBL.rev.1.bt2 ./index/GRCm38_101_ENSEMBL.rev.1.bt2
mv GRCm38_101_ENSEMBL.rev.2.bt2 ./index/GRCm38_101_ENSEMBL.rev.2.bt2

sleep 1

rm Mus_musculus.GRCm38.cdna.all.fa
rm Mus_musculus.GRCm38.ncrna.fa
rm Mouse_RNA.fa
rm Mus_musculus.GRCm38.rna.all.fa

