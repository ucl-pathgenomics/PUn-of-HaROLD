# PUN of HaROLD: Pretty Useful ruNtime of HaROLD

This program automates the running of processed in performing haplotype reconstruction on longitudinal deep sequencing samples, by analysing co-varying variants in a probabilistic framework. 

## Instructions:
-requires n bam files are stored in the /01_bam folder
-requires a .txt file in /01_bam that has /n teparated filenames (to define the order of processing)
-requires a .fasta file in /01_bam that is the reference file which the .bam files were assembled against
-define the commands as in the top instruciton line in the python file


requirments:
- packaged in here is a punofharold.yml file, use this to create a conda environment for runtime.
