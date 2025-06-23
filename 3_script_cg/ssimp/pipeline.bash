#!/bin/bash

## munge summary statistics into format for SSIMP
perl munge_for_ssimp.pl ../DATASETS/MR_DATA_V2/hcmr_NTproBNP_gwas_ssf.h.tsv.gz

## Run SSIMP
sbatch --array=2-224 run_ssimp_mr2.sh

## Liftover to build 38
sbatch --array=3-224:5 liftover_idx_mr.sh

## munge Imputed variants back into original summstats
perl munge_into_original_mr.pl tadros24_hcm_gwas_ssf.h.tsv.gz
