#!/bin/bash
#SBATCH -p long
#SBATCH -A procardis.prj
#SBATCH --constraint="skl-compat"
#SBATCH --mem=30G

# nextflow run EBISPOT/gwas-sumstats-harmoniser -r v1.1.10 \
#     --ref '/well/PROCARDIS/jchan/bin/gwas-sumstats-harmoniser/reference_files/' \
#     --file '/well/PROCARDIS/jchan/misc/gwas_summstats/HCM/Tadros24/nonMTAG/hcm.fix.2024.format_gwas_ssf.tsv.gz' \
#     --harm \
#     --chromlist '1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22' \
#     --profile executor,conda \
#     --to_build '38' 


nextflow run EBISPOT/gwas-sumstats-harmoniser -r v1.1.10 \
    --ref '/well/PROCARDIS/jchan/bin/gwas-sumstats-harmoniser/reference_files/' \
    --file '/well/PROCARDIS/jchan/misc/gwas_summstats/CMR/rvsvi/pirrucello22/pirrucello22_rvsvi_gwas_ssf.tsv.gz' \
    --harm \
    --chromlist '1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22' \
    --profile executor,conda \
    --to_build '38' 
