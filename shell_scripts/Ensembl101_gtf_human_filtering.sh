#!/bin/bash -i    
# Genome metadata
# oridinal : https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr
# https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#hg19_3.0.0
# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
# - These filters are relevant only to GTF files from GENCODE. The GTFs from
#   Ensembl release 98 have the following differences:
#   - The names "gene_biotype" and "transcript_biotype" are used instead of
#     "gene_type" and "transcript_type".
#   - Readthrough transcripts are present but are not marked with the
#     "readthrough_transcript" tag.

wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz
unpigz Homo_sapiens.GRCh38.101.gtf.gz

BIOTYPE_PATTERN=\
"(protein_coding|lincRNA|antisense|bidirectional_promoter_lnkRNA|macro_lnkRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_biotype \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_biotype \"${BIOTYPE_PATTERN}\""
#READTHROUGH_PATTERN="tag \"readthrough_transcript\""


# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
cat Homo_sapiens.GRCh38.101.gtf \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > Homo_sapiens.GRCh38.101.gene_allowlist.txt

#    | grep -Ev "$READTHROUGH_PATTERN" \

#split gene_allowlist file for parallel processing
split -1947 Homo_sapiens.GRCh38.101.gene_allowlist.txt Homo_sapiens.GRCh38.101.gene_allowlist-

# Filter the GTF file based on the gene allowlist
ls Homo_sapiens.GRCh38.101.gene_allowlist-* | cat | parallel -P 16 -a \
- 'sh ./shell_scripts/gtf_filtering.sh' {} Homo_sapiens.GRCh38.101.gtf 

ls Homo_sapiens.GRCh38.101.gene_allowlist-* | cat | xargs cat > Homo_sapiens.GRCh38.101.gene_filtered1.gtf
rm Homo_sapiens.GRCh38.101.gene_allowlist-*

# Copy header lines beginning with "#"
grep -E "^#" Homo_sapiens.GRCh38.101.gtf > Homo_sapiens.GRCh38.101.header.gtf
cat Homo_sapiens.GRCh38.101.header.gtf Homo_sapiens.GRCh38.101.gene_filtered1.gtf > Homo_sapiens.GRCh38.101.gene_filtered.gtf 
rm Homo_sapiens.GRCh38.101.header.gtf
rm Homo_sapiens.GRCh38.101.gene_filtered1.gtf 
rm Homo_sapiens.GRCh38.101.gene_allowlost.txt
pigz Homo_sapiens.GRCh38.101.gene_filtered.gtf
mv Homo_sapiens.GRCh38.101.gene_filtered.gtf.gz ./reference/fasta_reference_human/Homo_sapiens.GRCh38.101.gene_filtered.gtf.gz
rm Homo_sapiens.GRCh38.101.gtf

