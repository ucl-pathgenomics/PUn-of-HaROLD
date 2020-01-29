
#!/bin/bash -l
#$ -S /bin/bash
#$ -o C:\Oscar\OneDrive\UCL\19-5_punofharold\refinement\file-0/out
#$ -e C:\Oscar\OneDrive\UCL\19-5_punofharold\refinement\file-0/error
#$ -l h_rt=24:00:00
#$ -l tmem=12.9G,h_vmem=12.9G
#$ -N  h_6
#$ -wd C:\Oscar\OneDrive\UCL\19-5_punofharold\refinement\file-0/
#$ -V
#$ -R y



### activate conda
source /share/apps/anaconda/bin/activate /home/ocharles/.conda/envs/harold-oct19

### HAROLD
cd C:\Oscar\OneDrive\UCL\19-5_punofharold\refinement\file-0/
java -Xms1500m -Xmx11000m -jar C:\Oscar\OneDrive\UCL\19-5_punofharold\03_dependencies/HAROLD/Cluster_RG/dist/Cluster_RG.jar --count-file  --haplotypes 6 --alpha-frac 0.01 --gamma-cache 10000 -H -L --threads 0 -f C:\Oscar\OneDrive\UCL\19-5_punofharold/h5_a0.01/freq_table.tsv > C:\Oscar\OneDrive\UCL\19-5_punofharold\refinement\file-0/out/harold_log.txt

### end of script messages
echo "********************************* SCRIPT COMPLETED *********************************"

