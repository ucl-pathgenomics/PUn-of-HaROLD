#!/bin/bash -l
#$ -S /bin/bash
#$ -o ${dir_curr}/out
#$ -e ${dir_curr}/error
#$ -l h_rt=64:00:00
#$ -l tmem=12.9G,h_vmem=12.9G
#$ -N  refinement
#$ -wd ${dir_proj}/
#$ -V
#$ -R y

mkdir ${dir_curr}/out
mkdir ${dir_curr}/error


### activate conda
source ${conda_activate_loc} ${conda_env}

cd ${dir_proj}/

### Refinement - we just need the Bam from one of the files tested to refine, not a list of files
${v_java} -cp ${dir_harold}/HAROLD/Cluster_RG/dist/lib/htsjdk-unspecified-SNAPSHOT.jar:${dir_harold}/HAROLD/Cluster_RG/dist/lib/pal-1.5.1.jar:${dir_harold}/HAROLD/Cluster_RG/dist/lib/commons-math3-3.5.jar:${dir_harold}/HAROLD/RefineHaplotypes/dist/RefineHaplotypes.jar refineHaplotypes.RefineHaplotypes -tag output.txt -hapSeq ${refine.best_harold_dir}/ResultHaplo.fasta -refSeq ${ref.fasta} -baseFreq ${refine.best_harold_dir}/Results.lld -bam ${yoursample.bam} -iterate > ${dir_curr}/RefineResults.txt

##java -cp 03_dependencies/HAROLD/RefineHaplotypes/dist/RefineHaplotypes.jar:03_dependencies/HAROLD/Cluster_RG/dist/lib/htsjdk-unspecified-SNAPSHOT.jar:03_dependencies/HAROLD/Cluster_RG/dist/lib/pal-1.5.1.jar:03_dependencies/HAROLD/Cluster_RG/dist/lib/commons-math3-3.5.jar refineHaplotypes.RefineHaplotypes -tag refinement/output.txt -hapSeq h6_a0.01/ResultsHaplo.fasta -refSeq 01_bam/ref.fasta -baseFreq h6_a0.01/Results.lld -bam 01_bam/CMV-2064_nodups_sorted.bam -iterate

### end of script messages
echo "********************************* SCRIPT COMPLETED *********************************"