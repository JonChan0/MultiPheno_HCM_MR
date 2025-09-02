#!/bin/bash

## munge to coloc format
perl munge_withpos.pl thanaj22_radialdsr_gwas_ssf.h.with_n.tsv.gz quant radialdsr

## run coloc
sbatch --array=2-192:20 run_coloc_v2.sh

## generate coloc output
sbatch --array=2-192:10 parse_rds_v2.sh

## munge summary of results
perl summary.pl > coloc_nonMTAG.out
