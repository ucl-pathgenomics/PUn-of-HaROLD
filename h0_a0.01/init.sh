# pipeline for pun of harold using cluster
# a template for replacing items with varables

v_threadnum=`cat C:\Oscar\OneDrive\UCL\19-5_punofharold/h0_a0.01\pat_14_harold_in.txt | wc -l`
mkdir C:\Oscar\OneDrive\UCL\19-5_punofharold/out
mkdir C:\Oscar\OneDrive\UCL\19-5_punofharold/error

################## pipeline #################

    echo "
#!/bin/bash -l
#$ -S /bin/bash
#$ -o C:\Oscar\OneDrive\UCL\19-5_punofharold/out
#$ -e C:\Oscar\OneDrive\UCL\19-5_punofharold/error
#$ -l h_rt=24:00:00
#$ -l tmem=12.9G,h_vmem=12.9G
#$ -N  h_0
#$ -wd C:\Oscar\OneDrive\UCL\19-5_punofharold/
#$ -V
#$ -R y

### activate conda
source /share/apps/anaconda/bin/activate /home/ocharles/.conda/envs/harold-oct19
### run diversity?


### HAROLD
cd C:\Oscar\OneDrive\UCL\19-5_punofharold
java -jar C:\Oscar\OneDrive\UCL\19-5_punofharold/03_dependencies/HAROLD/Cluster_RG/dist/Cluster_RG.jar --count-file C:\Oscar\OneDrive\UCL\19-5_punofharold/h0_a0.01\pat_14_harold_in.txt --haplotypes 0 --alpha-frac 0.01 --gamma-cache 10000 -H -L --threads 6  > out/harold_log.txt

### end of script messages
echo "********************************* SCRIPT COMPLETED *********************************"


 " > C:\Oscar\OneDrive\UCL\19-5_punofharold/hap_0.sh

     qsub C:\Oscar\OneDrive\UCL\19-5_punofharold/hap_0.sh

