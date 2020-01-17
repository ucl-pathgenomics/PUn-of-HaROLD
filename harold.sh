    # pipeline for pun of harold using cluster
# a template for replacing items with varables

mkdir ${dir_curr}/out
mkdir ${dir_curr}/error

################## pipeline #################

    echo "
#!/bin/bash -l
#$ -S /bin/bash
#$ -o ${dir_curr}/out
#$ -e ${dir_curr}/error
#$ -l h_rt=24:00:00
#$ -l tmem=12.9G,h_vmem=12.9G
#$ -N  h_${v_haplonum}
#$ -wd ${dir_curr}/
#$ -V
#$ -R y

### activate conda
source ${conda_activate_loc} ${conda_env}

### HAROLD
cd ${dir_curr}/
${v_java} -jar ${dir_harold}/HAROLD/Cluster_RG/dist/Cluster_RG.jar --count-file ${infile} --haplotypes ${v_haplonum} --alpha-frac ${v_alphafrac} --gamma-cache 10000 -H -L --threads ${v_threadnum} ${init_freq} > ${dir_curr}/out/harold_log.txt

### end of script messages
echo "********************************* SCRIPT COMPLETED *********************************"


 " > ${dir_curr}/init.sh

     ${command} ${dir_curr}/init.sh

