#release-101 build
#build date 20210127
#must set hisat2 path

wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/snp142Common.txt.gz
wget ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
sleep 1

mv Mus_musculus.GRCm38.dna.primary_assembly.fa.gz Mus_musculus.GRCm38.101.dna.primary_assembly.fa.gz
cp Mus_musculus.GRCm38.101.dna.primary_assembly.fa.gz ./reference/fasta_reference_mouse/Mus_musculus.GRCm38.101.dna.primary_assembly.fa.gz
cp ./reference/fasta_reference_mouse/Mus_musculus.GRCm38.101.gene_filtered.gtf.gz Mus_musculus.GRCm38.101.gene_filtered.gtf.gz
unpigz Mus_musculus.GRCm38.101.gene_filtered.gtf.gz 
unpigz Mus_musculus.GRCm38.101.dna.primary_assembly.fa.gz
unpigz snp142Common.txt.gz
sleep 1

awk 'BEGIN{OFS="\t"} {if($2 ~ /^chr/) {$2 = substr($2, 4)}; if($2 == "M") {$2 = "MT"} print}' snp142Common.txt > snp142Common.txt.ensembl
hisat2_extract_snps_haplotypes_UCSC.py Mus_musculus.GRCm38.dna.primary_assembly.fa snp142Common.txt.ensembl genome
hisat2_extract_splice_sites.py Mus_musculus.GRCm38.101.gene_filtered.gtf > genome.ss
hisat2_extract_exons.py Mus_musculus.GRCm38.101.gene_filtered.gtf > genome.exon


threads=`grep -ce '^processor\s\+:' /proc/cpuinfo`
tmp=`echo $((threads/2))`
hisat2-build -p ${tmp} Mus_musculus.GRCm38.101.dna.primary_assembly.fa --snp genome.snp --haplotype genome.haplotype --ss genome.ss --exon genome.exon GRCm38_101_ENSEMBL

sleep 1

mv GRCm38_101_ENSEMBL.1.ht2 ./index/GRCm38_101_ENSEMBL.1.ht2
mv GRCm38_101_ENSEMBL.2.ht2 ./index/GRCm38_101_ENSEMBL.2.ht2
mv GRCm38_101_ENSEMBL.3.ht2 ./index/GRCm38_101_ENSEMBL.3.ht2
mv GRCm38_101_ENSEMBL.4.ht2 ./index/GRCm38_101_ENSEMBL.4.ht2
mv GRCm38_101_ENSEMBL.5.ht2 ./index/GRCm38_101_ENSEMBL.5.ht2
mv GRCm38_101_ENSEMBL.6.ht2 ./index/GRCm38_101_ENSEMBL.6.ht2
mv GRCm38_101_ENSEMBL.7.ht2 ./index/GRCm38_101_ENSEMBL.7.ht2
mv GRCm38_101_ENSEMBL.8.ht2 ./index/GRCm38_101_ENSEMBL.8.ht2

rm genome.ss
rm genome.haplotype
rm genome.snp
rm genome.exon
rm snp142Common.txt.ensembl
rm snp142Common.txt
rm Mus_musculus.GRCm38.101.dna.primary_assembly.fa
rm Mus_musculus.GRCm38.101.gene_filtered.gtf

