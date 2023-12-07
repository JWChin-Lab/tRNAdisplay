# tRNA display analysis

Scripts associated for analysis of Illumina NGS data from tRNA display experiments, associated with Dunkelmann, Piedrafita et al., Nature 2023).  

Illumina paired end fastq.gz files are first merged with PEAR (Zhang et al., Bioinformatics 2014) and aligned with BowTie2 (Langmead et al., Nature Methods 2012).

Resulting sam files (as .txt files) are translated and processed using custom R scripts:
\ 1_  
\ 2_
\ 3_ 

The two raw paired-end read files per sample were merged using PEAR (Zhang et al. 2014) using the command:

the resulting text file was aligned to a reference sequence (ref.a) using Bowtie2 (Langmead et al. 2012) using the command:

The resulting text file was translated and processed using script XXX to yield a list of residue combinations at the targeted positions.

These frequency lists were used as input for the tRNA display analysis script (script YYY) yielding spindle plots and gated hit lists.

