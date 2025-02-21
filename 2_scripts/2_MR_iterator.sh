#This script iterates over all the config files in ../input/configs to run the mendelian_randomisation.smk pipeline for each pair of traits

module load Anaconda3/2022.05
eval "$(conda shell.bash hook)"
conda activate gms

echo $(date)


#This applies to running the bidirectional MR with HCMR GWAS of plasma proteins with HCM
for file in ../1_data/input_configs/*.yaml
do
    echo "Carrying out bidirectional MR on ${file}"
    snakemake --profile bmrc_profile_smk5  --snakefile 2_bidirectional_MR.smk --configfile $file
done








