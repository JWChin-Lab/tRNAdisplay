# tRNA display analysis

Scripts associated for analysis of Illumina NGS data from tRNA display experiments, associated with Dunkelmann, Piedrafita et al., Nature 2023).  

Illumina paired end fastq.gz files are first merged with PEAR (Zhang et al., Bioinformatics 2014) and aligned with BowTie2 (Langmead et al., Nature Methods 2012):



Resulting sam files (as .txt files) are translated and processed using custom R scripts:  \
1_extract_and_sum_variants.R filters and counts sequences.\
2_translate_variants.R translates and yields residue combinations at targeted positions. \
3_tRNA_display_analysis.R plots spindle plots and produces gated hit lists.  


