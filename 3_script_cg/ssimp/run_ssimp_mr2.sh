#!/bin/bash 

#SBATCH -p short
#SBATCH -o slurm/ssimp_%A.%a.out
#SBATCH -e slurm/ssimp_%A.%a.err
 
echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "------------------------------------------------" 
 
date

FIRST=$SLURM_ARRAY_TASK_ID

LAST=$(($FIRST + $SLURM_ARRAY_TASK_STEP - 1))

for (( MY_TASK=$FIRST; MY_TASK<=$LAST; MY_TASK++ ))
do
    chr=`grep -w ^$MY_TASK nonMTAG_34loci_indexes.txt | cut -f 2`
    st=`grep -w ^$MY_TASK nonMTAG_34loci_indexes.txt | cut -f 3`
    en=`grep -w ^$MY_TASK nonMTAG_34loci_indexes.txt | cut -f 4`
    file=`grep -w ^$MY_TASK nonMTAG_34loci_indexes.txt | cut -f 5`

    bash run_ssimp.bash $chr $st $en $file
done

date
