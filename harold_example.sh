# have confirmed for CMV alpha frac of 0.01 is the same performance as 0.1
#pipeline variables
project=19_4-CMV
patient=17
v_haplonum=3
v_alphafrac=0.01
dir_harold=/cluster/scratch8b/oscar/RichardProgram
dir_out=/cluster/scratch8b/oscar/${project}/pat_${patient}
infile1=${dir_out}/pat_${patient}_harold_in.txt
infile2=${dir_out}/haplo_${v_haplonum}_alpha_${v_alphafrac}/pat_${patient}_harold_in.txt
mkdir ${dir_out}
mkdir ${dir_out}/haplo_${v_haplonum}_alpha_${v_alphafrac}
mkdir ${dir_out}/haplo_${v_haplonum}_alpha_${v_alphafrac}/out
mkdir ${dir_out}/haplo_${v_haplonum}_alpha_${v_alphafrac}/error


################## setup ###################
#copy strancount files
echo " /cluster/scratch8b/oscar/19_4-CMV/CMV-2096/refbased/CMV-2096_nodups_sorted_diversity.strandcount.csv
/cluster/scratch8b/oscar/19_4-CMV/CMV-2098/refbased/CMV-2098_nodups_sorted_diversity.strandcount.csv
/cluster/scratch8b/oscar/19_4-CMV/CMV-2103/refbased/CMV-2103_nodups_sorted_diversity.strandcount.csv
/cluster/scratch8b/oscar/19_4-CMV/CMV-2104/refbased/CMV-2104_nodups_sorted_diversity.strandcount.csv
/cluster/scratch8b/oscar/19_4-CMV/CMV-2108/refbased/CMV-2108_nodups_sorted_diversity.strandcount.csv
/cluster/scratch8b/oscar/19_4-CMV/CMV-2145/refbased/CMV-2145_nodups_sorted_diversity.strandcount.csv" > ${infile1}


xargs -a ${infile1} cp -t ${dir_out}/haplo_${v_haplonum}_alpha_${v_alphafrac}/
#create and harold required file list
echo "CMV-2096_nodups_sorted_diversity.strandcount.csv
CMV-2098_nodups_sorted_diversity.strandcount.csv
CMV-2103_nodups_sorted_diversity.strandcount.csv
CMV-2104_nodups_sorted_diversity.strandcount.csv
CMV-2108_nodups_sorted_diversity.strandcount.csv
CMV-2145_nodups_sorted_diversity.strandcount.csv" > ${infile2}

v_threadnum=`cat ${dir_out}/haplo_${v_haplonum}_alpha_${v_alphafrac}/pat_${patient}_harold_in.txt | wc -l`

################## pipeline #################

    echo "
#!/bin/bash -l
#$ -S /bin/bash
#$ -o ${dir_out}/haplo_${v_haplonum}_alpha_${v_alphafrac}/out
#$ -e ${dir_out}/haplo_${v_haplonum}_alpha_${v_alphafrac}/error
#$ -l h_rt=24:00:00
#$ -l tmem=12.9G,h_vmem=12.9G
#$ -N  h_${patient}_${v_haplonum}
#$ -wd ${dir_out}/haplo_${v_haplonum}_alpha_${v_alphafrac}/
#$ -V
#$ -R y


## activate conda
source /share/apps/anaconda/bin/activate /home/ocharles/.conda/envs/harold-oct19

#list input files
# find ${dir_in}/${project}/ -name "*strandcount*" > all_strandcount.txt

### HAROLD
cd ${dir_out}/haplo_${v_haplonum}_alpha_${v_alphafrac}/
java -jar ${dir_harold}/HAROLD/Cluster_RG/dist/Cluster_RG.jar --count-file ${infile2} --haplotypes ${v_haplonum} --alpha-frac ${v_alphafrac} --gamma-cache 10000 -H -L --threads ${v_threadnum} > harold_log.txt

#Note: -threads = parallel run, set this as the number of samples you have in each run;
#-haplotypes = expected number of haplotypes 
#When you rerun HAROLD using other expected number of haplotypes, try to specify the frequencies in a table and feed it to the program using --initial-freq-file


 
### refinement
#java -Xmx1000m -cp ${dir_harold}/Documents/HAROLD/Cluster_RG/dist/lib/htsjdk-unspecified-SNAPSHOT.jar:${dir_harold}/HAROLD/Cluster_RG/dist/lib/pal-1.5.1.jar:${dir_harold}/HAROLD/Cluster_RG/dist/lib/commons-math3-3.5.jar:${dir_harold}/HAROLD/RefineHaplotypes/dist/RefineHaplotypes.jar refineHaplotypes.RefineHaplotypes -tag nameofoutputfile -hapSeq outputseqfromHAROLD -refSeq ref.fasta -baseFreq lldfilefromHAROLD -bam yoursample.bam -iterate > RefineResults_yoursample





### end of script messages
echo "********************************* SCRIPT COMPLETED *********************************"
date
echo "********************************* SCRIPT COMPLETED *********************************"


 " > ${dir_out}/haplo_${v_haplonum}_alpha_${v_alphafrac}/${patient}_${v_haplonum}hap.sh

     qsub ${dir_out}/haplo_${v_haplonum}_alpha_${v_alphafrac}/${patient}_${v_haplonum}hap.sh

