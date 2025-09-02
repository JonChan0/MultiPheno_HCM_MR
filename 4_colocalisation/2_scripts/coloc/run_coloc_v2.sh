#!/bin/bash 

#SBATCH -p short
#SBATCH -o slurm/munge_%A.%a.out
#SBATCH -e slurm/munge_%A.%a.err
 
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
    loci=`grep -w ^$MY_TASK indexes.txt | cut -f 2`
    t1=`grep -w ^$MY_TASK indexes.txt | cut -f 3`
    t2=`grep -w ^$MY_TASK indexes.txt | cut -f 4`
    Rscript run_coloc_v2.R $loci $t1 $t2
done

date
