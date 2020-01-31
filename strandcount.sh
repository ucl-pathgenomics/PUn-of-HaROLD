#!/bin/bash -l
#$ -S /bin/bash
#$ -o ${dir_proj}/01_bam
#$ -e ${dir_proj}/01_bam
#$ -l h_rt=3:00:00
#$ -l tmem=20G,h_vmem=20G
#$ -N  bam_proc
#$ -wd ${dir_proj}/
#$ -V
#$ -R y

date
echo "********************************* SCRIPT STARTED *********************************"

### activate conda
source ${conda_activate_loc} ${conda_env}

cd ${dir_proj}/

### loop over each bam file

for file in 01_bam/*.bam; do
    echo file - ${file}
    out1=${file%.*}
    #echo ${out1}

    samtools view -h -b -G69 ${file} | samtools view -h -b -G133 > ${out1}_ready.bam
    
    java -cp 03_dependencies/HAROLD/Cluster_RG/dist/lib/htsjdk-unspecified-SNAPSHOT.jar:03_dependencies/HAROLD/Cluster_RG/dist/lib/picocli-3.6.0.jar:03_dependencies/HAROLD/Cluster_RG/dist/lib/pal-1.5.1.jar:03_dependencies/HAROLD/Cluster_RG/dist/lib/cache2k-all-1.0.2.Final.jar:03_dependencies/HAROLD/Cluster_RG/dist/lib/commons-math3-3.5.jar:03_dependencies/MakeReadCount.jar makereadcount.MakeReadCount ${out1}_ready.bam
    #move files once created
    mv *.log 02_strandcount
    mv *.csv 02_strandcount

    #delete intermediate bam
    rm 01_bam/*_ready.bam
done


### end of script messages
echo "********************************* SCRIPT COMPLETED *********************************"
date