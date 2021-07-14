#release-101 build
#build date 20210127
#must set hisat2 path

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/snp144Common.txt.gz
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
sleep 1

cp Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ./reference/fasta_reference_human/Homo_sapiens.GRCh38.101.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz Homo_sapiens.GRCh38.101.dna.primary_assembly.fa.gz
cp ./reference/fasta_reference_human/Homo_sapiens.GRCh38.101.gene_filtered.gtf.gz Homo_sapiens.GRCh38.101.gene_filtered.gtf.gz
unpigz Homo_sapiens.GRCh38.101.gene_filtered.gtf.gz
unpigz Homo_sapiens.GRCh38.101.dna.primary_assembly.fa.gz
unpigz snp144Common.txt.gz
sleep 1

awk 'BEGIN{OFS="\t"} {if($2 ~ /^chr/) {$2 = substr($2, 4)}; if($2 == "M") {$2 = "MT"} print}' snp144Common.txt > snp144Common.txt.ensembl
hisat2_extract_snps_haplotypes_UCSC.py Homo_sapiens.GRCh38.101.dna.primary_assembly.fa snp144Common.txt.ensembl genome
hisat2_extract_splice_sites.py Homo_sapiens.GRCh38.101.gene_filtered.gtf > genome.ss
hisat2_extract_exons.py Homo_sapiens.GRCh38.101.gene_filtered.gtf > genome.exon


threads=`grep -ce '^processor\s\+:' /proc/cpuinfo`
tmp=`echo $((threads/2))`

hisat2-build -p ${tmp} Homo_sapiens.GRCh38.101.dna.primary_assembly.fa --snp genome.snp --haplotype genome.haplotype --ss genome.ss --exon genome.exon GRCh38_101_ENSEMBL

mv GRCh38_101_ENSEMBL.1.ht2 ./index/GRCh38_101_ENSEMBL.1.ht2
mv GRCh38_101_ENSEMBL.2.ht2 ./index/GRCh38_101_ENSEMBL.2.ht2
mv GRCh38_101_ENSEMBL.3.ht2 ./index/GRCh38_101_ENSEMBL.3.ht2
mv GRCh38_101_ENSEMBL.4.ht2 ./index/GRCh38_101_ENSEMBL.4.ht2
mv GRCh38_101_ENSEMBL.5.ht2 ./index/GRCh38_101_ENSEMBL.5.ht2
mv GRCh38_101_ENSEMBL.6.ht2 ./index/GRCh38_101_ENSEMBL.6.ht2
mv GRCh38_101_ENSEMBL.7.ht2 ./index/GRCh38_101_ENSEMBL.7.ht2
mv GRCh38_101_ENSEMBL.8.ht2 ./index/GRCh38_101_ENSEMBL.8.ht2

rm genome.ss
rm genome.haplotype
rm genome.snp
rm genome.exon
rm snp144Common.txt.ensembl
rm snp144Common.txt
rm Homo_sapiens.GRCh38.101.dna.primary_assembly.fa
rm Homo_sapiens.GRCh38.101.gene_filtered.gtf

