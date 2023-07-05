#!/bin/bash

# Load the required modules
module load bwa/0.7.17
module load samtools/1.14
module load bcftools/1.15

# Create a working directory
mkdir -p /project/pi_frederic_chain_uml_edu/Hackbio/dc_workshop
cd /project/pi_frederic_chain_uml_edu/Hackbio/dc_workshop

# Download the stickleback reference genome
curl -L -o stickleback.fa.gz https://stickleback.genetics.uga.edu/downloadData/v5.0.1_assembly/stickleback_v5.0.1_assembly.fa.gz

# Uncompress the reference genome file
gunzip stickleback.fa.gz

# Create a directory for the results
mkdir -p results

# Download the FastQ files for the sample
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR770/ERR770591/ERR770591_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR770/ERR770591/ERR770591_2.fastq.gz

# Generate index files
bwa index stickleback.fa

# Align reads to the reference genome using BWA
bwa mem stickleback.fa ERR770591_1.fastq.gz ERR770591_2.fastq.gz > results/aligned.sam

# Convert the SAM file to BAM format using samtools
samtools view -S -b results/aligned.sam > results/aligned.bam

# Sort the BAM file
samtools sort results/aligned.bam -o results/sorted.bam

# Calculate the read coverage of positions in the genome
bcftools mpileup -O b -o results/variants.bcf -f stickleback.fa --threads 8 -q 20 -Q 30 results/sorted.bam

# Detect the single nucleotide variants (SNVs)
bcftools call --ploidy 1 -m -v -o results/variants.vcf results/variants.bcf

# Filter and report the SNV variants in variant calling format (VCF)
bcftools view results/variants.vcf -o results/filtered_variants.vcf

# Explore the VCF format
less -S results/filtered_variants.vcf

# Use the grep and wc commands to assess how many variants are in the VCF file
grep -v "#" results/filtered_variants.vcf | wc -l

# Optional step: Assess the alignment (visualization)
samtools index results/sorted.bam

# Viewing with tview
samtools tview results/sorted.bam stickleback.fa

# Using the 'cat' command to concatenate the lines of the variant view and then pipe the output to a text file named 'variant_view.txt'
cat << EOF > variant_view.txt
1         11        21        31        41        51        61        71        81        91        101       111       121       131       141       151
CTTCTTCTCCTCTTCCTCTTCCTTCTCCCTCTTCCTCCCTGTCGGCGTGTCATCAGATCTGACCAGTGTGTGTGTGTGTTTGCGTGCGTGTGCGTGCGTGCGTGCGTGTTTGCGTGCGTGTTTGCGTGCGTGTTTGCGTGCGTGTGTGTGTGTGTGTGTGTGTG
EOF

# To view the content of 'variant_view.txt'
cat variant_view.txt