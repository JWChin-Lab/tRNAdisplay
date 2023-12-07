# tRNA display analysis

Scripts associated for analysis of Illumina NGS data from tRNA display experiments, associated with Dunkelmann, Piedrafita et al., Nature 2023).  

Illumina paired end fastq.gz files are first merged with PEAR (Zhang et al., Bioinformatics 2014) and aligned with BowTie2 (Langmead et al., Nature Methods 2012). 

#!/bin/bash \
(pear -f $READS1 -r $READS2 -o 1_joined) > 1_PEAR_summary.txt \
bowtie2-build ref1.fa ref1 \
(bowtie2 -x ref1 -U 1_joined.txt.assembled.fastq -S 2_aligned.txt) 2>2_alignmentsummary.txt 

Resulting sam files (as .txt files) are translated and processed using custom R scripts: 

1_extract_and_sum_variants.R filters and counts sequences.\
2_translate_variants.R translates and yields residue combinations at targeted positions. \
3_tRNA_display_analysis.R plots spindle plots and produces gated hit lists. Requires testthat, ggplot2, overlapping and factoextra



