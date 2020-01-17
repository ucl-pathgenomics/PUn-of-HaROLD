
#!/bin/bash -l
#$ -S /bin/bash
#$ -o /cluster/scratch8b/oscar/19_4-CMV/pat_14/haplo_6_alpha_0.01/out
#$ -e /cluster/scratch8b/oscar/19_4-CMV/pat_14/haplo_6_alpha_0.01/error
#$ -l h_rt=24:00:00
#$ -l tmem=12.9G,h_vmem=12.9G
#$ -N  harold_14_6
#$ -wd /cluster/scratch8b/oscar/19_4-CMV/pat_14/haplo_6_alpha_0.01/
#$ -V
#$ -R y


## activate conda
source /share/apps/anaconda/bin/activate /home/ocharles/.conda/envs/harold-oct19

#list input files
# find /cluster/scratch8b/oscar/19_4-CMV/ -name *strandcount* > all_strandcount.txt

### HAROLD
cd /cluster/scratch8b/oscar/19_4-CMV/pat_14/haplo_6_alpha_0.01/
java -jar /cluster/scratch8b/oscar/RichardProgram/HAROLD/Cluster_RG/dist/Cluster_RG.jar --count-file /cluster/scratch8b/oscar/19_4-CMV/pat_14/pat_14_harold_in.txt --haplotypes 6 --alpha-frac 0.01 --gamma-cache 10000 -H -L --threads 5 > harold_log.txt

#Note: -threads = parallel run, set this as the number of samples you have in each run;
#-haplotypes = expected number of haplotypes 
#When you rerun HAROLD using other expected number of haplotypes, try to specify the frequencies in a table and feed it to the program using --initial-freq-file


 
### refinement
#java -Xmx1000m -cp /cluster/scratch8b/oscar/RichardProgram/Documents/HAROLD/Cluster_RG/dist/lib/htsjdk-unspecified-SNAPSHOT.jar:/cluster/scratch8b/oscar/RichardProgram/HAROLD/Cluster_RG/dist/lib/pal-1.5.1.jar:/cluster/scratch8b/oscar/RichardProgram/HAROLD/Cluster_RG/dist/lib/commons-math3-3.5.jar:/cluster/scratch8b/oscar/RichardProgram/HAROLD/RefineHaplotypes/dist/RefineHaplotypes.jar refineHaplotypes.RefineHaplotypes -tag nameofoutputfile -hapSeq outputseqfromHAROLD -refSeq ref.fasta -baseFreq lldfilefromHAROLD -bam yoursample.bam -iterate > RefineResults_yoursample





### end of script messages
echo ********************************* SCRIPT COMPLETED *********************************
date
echo ********************************* SCRIPT COMPLETED *********************************


 
