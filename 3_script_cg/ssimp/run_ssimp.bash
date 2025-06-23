#!/bin/bash

i=$1
s=$2
e=$3
file=$4

#echo $i $s $e

root="/well/PROCARDIS/cgrace/ssimp"

${root}/ssimp --gwas ${file} \
    --ref ${root}/1000g_p3_all/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
    --out chunks_mr_nonMTAG_34loci/${file}.chr${i}.${s}.${e} \
    --impute.range ${i}:${s}-${i}:${e}  \
    --log chunks_mr_nonMTAG_34loci/${file}.chr${i}.${s}.${e}.log
