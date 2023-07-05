# Chain's Lab Variant Calling

Title: Laboratory Report - Variant Calling Analysis of Stickleback Fish Genome
Abstract: This laboratory report presents the experimental workflow and results of a variant calling analysis conducted on the genome data of three-spined stickleback fish (Gasterosteus aculeatus). The study aimed to identify single nucleotide variants (SNVs) in the stickleback genome using high-throughput sequencing data. The analysis pipeline involved read alignment, BAM file sorting, variant calling, and filtering steps. The work was performed on the Unity Cluster at the University of Massachusetts Lowell under the supervision of Dr. Frédéric JJ Chain.
Introduction: The three-spined stickleback fish serves as a model organism for studying evolutionary genetics and adaptation. This laboratory study focused on analyzing the genome data of stickleback fish to identify genetic variations in the form of SNVs. The variant calling analysis aimed to provide insights into the genomic diversity and potential functional implications of these variants.
Methods:
Data Retrieval:
Two paired-end FastQ files (ERR770591_1.fastq.gz and ERR770591_2.fastq.gz) were obtained from the European Bioinformatics Institute (EBI) FTP server. These files corresponded to a transcriptome profiling study conducted by Huang et al. (2016) on stickleback fish immune tissues.
Preprocessing and Alignment:
The BWA software (version 0.7.17) was utilized to generate index files for the reference genome.
The FastQ files were aligned to the reference genome (stickleback.fa) using the BWA-MEM algorithm, resulting in a SAM file (aligned.sam).
File Conversion and Sorting:
The SAM file was converted to the compressed BAM format using SAMtools.
The BAM file (aligned.bam) was sorted by genomic coordinates using the SAMtools sort command.
Variant Calling and Filtering:
Variant calling was performed using BCFTOOLS mpileup to calculate read coverage and detect SNVs in the aligned reads. The results were saved in a BCF file (variants.bcf).
BCFTOOLS call was applied to the BCF file with a ploidy of 1 to generate a VCF file (variants.vcf) containing the identified SNVs.
The VCF file was further processed using BCFTOOLS view to filter and report the SNV variants in a separate VCF file (filtered_variants.vcf).

Results: The variant calling analysis identified a set of SNVs in the stickleback fish genome. The filtered_variants.vcf file contained the filtered and reported SNV variants. The analysis provided valuable insights into the genetic diversity and potential functional implications of these variants in stickleback fish.
Conclusion: This laboratory study successfully conducted a variant calling analysis on the genome data of three-spined stickleback fish. The analysis pipeline involved read alignment, BAM file sorting, variant calling, and filtering steps. The identification of SNV variants contributes to our understanding of the genetic variation present in stickleback fish and provides a foundation for further functional analysis. The findings pave the way for future investigations into the evolutionary genetics and adaptation of this model organism.
Reference: Huang, Y., Chain, F. J., Panchal, M., Eizaguirre, C., Kalbe, M., Lenz, T. L., ... & Feulner, P. G. (2016). Transcriptome profiling of immune tissues reveals habitat‐specific gene expression between lake and river sticklebacks. Molecular ecology, 25(4), 943-958.
