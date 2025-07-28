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
    bash munge_mr.bash $MY_TASK    
done

date
