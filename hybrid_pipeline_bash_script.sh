# Bash Script for Hybrid Human CoV2 Pipeline
# Created by Prof. Steve Shen and organized by Hezkiel Nanda
# Last update: 28-Aug-2020

# 1. Download and uncompress CoV2 genome and annotation files
# downloading genome sequencing and annotation (gff3) from ncbi
# https://www.ncbi.nlm.nih.gov/nuccore/MN985325.1?report=gbwithparts&log$=seqview

# for comparison: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512

# using cufflinks utility of MSI cluster to convert gff3 to gtf

# download MN985325.1 with Primary Assembly (GCA_009937904.1)
# from https://www.ncbi.nlm.nih.gov/assembly/GCA_009937905.1/

# fasta nucleodic acid
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/937/905/GCA_009937905.1_ASM993790v1/GCA_009937905.1_ASM993790v1_cds_from_genomic.fna.gz
# gff
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/937/905/GCA_009937905.1_ASM993790v1/GCA_009937905.1_ASM993790v1_genomic.gff.gz

# uncompress files
zcat GCA_009937905.1_ASM993790v1_cds_from_genomic.fna.gz > MN985325.fa
zcat GCA_009937905.1_ASM993790v1_genomic.gff.gz > MN985325.gff


# 2. Convert gff with Cufflinks (v2.2.1)
# http://cole-trapnell-lab.github.io/cufflinks/install/

gffread -E -T -O MN985325.gff -o MN985325.gtf

# -E  expose (warn about) duplicate transcript IDs and other potential
#      problems with the given GFF/GTF records
# -T  -o option will output GTF format instead of GFF3
# -O  process also non-transcript GFF records (by default non-transcript
#      records are ignored)
#  -o  the "filtered" GFF records will be written to <outfile.gff>
#      (use -o- for printing to stdout)

# Note: gtf might need to be modified to have an easier gene marker name


# 3. Download and uncompress Human genome and annotation
# Download all chromosome fasta/gtf files and merge/index to remove the unanchored scaffolds

# fasta
wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz 
# gtf
wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.chr.gtf.gz

# Note: we can use the latest release 

# merge and uncompress genome files
zcat *.fa.gz > Homo_sapiens.GRCh38.dna.chr.fa
zcat Homo_sapiens.GRCh38.98.chr.gtf.gz > Homo_sapiens.GRCh38.98.chr.gtf

# creating pre-mrna gtf file
# using formula to filter gtf file
awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript"{ $3="exon"; print}' Homo_sapiens.GRCh38.98.chr.gtf > Homo_sapiens.GRCh38.98.chr.premrna.gtf


# 4. Clean up gtf file with Cell Ranger (v3.0.1)
#https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

cellranger mkgtf Homo_sapiens.GRCh38.98.chr.premrna.gtf Homo_sapiens.GRCh38.98.chr.premrna_filt.gtf \
	--attribute=gene_biotype:protein_coding \
	--attribute=gene_biotype:lincRNA \
	--attribute=gene_biotype:antisense

# 5. Concat genome and annotation files
cat Homo_sapiens.GRCh38.dna.chr.fa MN985325.fa > GRCh38_cov2.fa

cat Homo_sapiens.GRCh38.98.chr.premrna_filt.gtf MN985325.gtf > GRCh38.98_cov2.gtf 


# 6. Create cell ranger genome reference
# Note: this requires high resources (ram and processing core)
cellranger mkref \
--genome=human_cov2_MN985325 \
--fasta=path/to/GRCh38_cov2.fa \
--genes=path/to/GRCh38.98_cov2.gtf


# 7. Process the sample data with the hybrid genome reference
# Note: this requires high resources (ram and processing core)
cellranger count \
	--id=project_id \
	--transcriptome=path/to/reference/human_cov2_MN985325 \
	--fastqs=path/to/fastqsample/ \
	--sample=sample_name \
	--expect-cells=1000