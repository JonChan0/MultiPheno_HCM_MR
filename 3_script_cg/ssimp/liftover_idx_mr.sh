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
#    chr=`grep -w ^$MY_TASK ../hcm_coord_b37.txt | cut -f 3`
#    st=`grep -w ^$MY_TASK ../hcm_coord_b37.txt | cut -f 4`
#    en=`grep -w ^$MY_TASK ../hcm_coord_b37.txt | cut -f 5`

#    bash run_ssimp.bash $chr $st $en sidorenko24_bmi_gwas_ssf.h.ssimp.tsv

    bash liftover_idx_mr_nonmtag.bash ${MY_TASK}
done

date
