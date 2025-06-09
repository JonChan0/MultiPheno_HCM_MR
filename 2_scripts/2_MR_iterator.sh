
#This script iterates over all the config files in ../input/configs to run the mendelian_randomisation.smk pipeline for each pair of traits

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate gms

echo $(date)

type=$1

#This applies to running the bidirectional MR with HCMR GWAS of plasma proteins with HCM
if [[ "$type" == 'bidirectional' ]]; then
    for file in ../1_data/input_configs/*.yaml
    do
        echo "Carrying out bidirectional MR on ${file}"
        snakemake --profile bmrc_profile_smk5  --snakefile 2_bidirectional_MR.smk --configfile $file --rerun-incomplete
    done

elif [[ "$type" == 'subgroup' ]]; then
    for file in ../1_data/input_configs/subgroup_analyses/*.yaml
    do
        echo "Carrying out UNIdirectional MR on ${file}"
        snakemake --profile bmrc_profile_smk5  --snakefile 2_unidirectional_MR.smk --configfile $file --rerun-incomplete
    done

elif [[ "$type" == 'sensitivity' ]]; then
    for file in ../1_data/input_configs/sensitivity/*.yaml
    do
        echo "Carrying out UNIdirectional MR on ${file}"
        snakemake --profile bmrc_profile_smk5  --snakefile 2_unidirectional_MR.smk --configfile $file --rerun-incomplete
    done

else
    echo "Please specify the type of MR you want to run: bidirectional, subgroup, or sensitivity"
    exit 1

fi

