#This script iterates over all the config files in ../input/configs to run the mendelian_randomisation.smk pipeline for each pair of traits

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate gms

echo $(date)


#This applies to running the bidirectional MR with HCMR GWAS of plasma proteins with HCM
for file in ../1_data/input_configs/subgroup_analyses/*.yaml
do
    echo "Carrying out UNI-directional MR on ${file}"
    snakemake --profile bmrc_profile_smk5  --snakefile 2_bidrectional_MR.smk --configfile $file --rerun-incomplete
done








