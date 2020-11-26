
# PUN of HaROLD: Pretty Useful ruNtime of HaROLD

NB: this is deprecated code, but some may find use of it.
  

This program automates the running of processed in performing haplotype reconstruction on longitudinal deep sequencing samples, by analysing co-varying variants in a probabilistic framework.

  

## Instructions:

 - Git clone this directory into folder
 - Save n bam files into <project_folder>/01_bam/
 - Save a **.txt** file including a newline separated list of your .bam files, in order of time (if you have it) into 01_bam/
 - Save a reference **.fasta**  file in 01_bam/
 - Edit pun_of_harold.py instructions:
	 - max-haplo - The maximum number of haplotypes to attempt in Haplotype_Reconstruction
	 - alpha-frac - The fraction of variable sites to include in error estimates,
		 -  0.01 is recommended. 0.1 for small genomes
		 - Larger is better, Larger is slower
	 - Command - the comm and requires to submit jobs in a cluster i.e. sh, qsub
	 - conda_activate_loc - the location of the binary to run anaconda
	 - conda_env - location of the conda environment to use in job submission
		 - punofharold.yml is compatible
	 - java - the command required to run java on your system


  
  
## requirements:

- packaged in here is a punofharold.yml file, use this to create a conda environment for runtime &/or job submission

## process:

 The process will run through from num_haplo =  2-max-haplo, running HAROLD.
 Formatiing outputs &  piping into the num_haplo + 1 run.
 Between num_haplo = 4 to max-haplo ifthe gradient of Log Likelihood for the most recent 3 runs is at any point negative; 
 HAROLD jobs will cease.
 A summary table & plot will be produced.
 For each BAM file in 01_bam/ haplotype refinement will be initialised, running each in parallel.
