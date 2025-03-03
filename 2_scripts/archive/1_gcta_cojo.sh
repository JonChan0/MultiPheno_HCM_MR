#!/bin/bash

sumstat_filepath = $1
output_basename = $2

for i in {1..22}
do
    /well/PROCARDIS/jchan/bin/gcta-1.94.3-linux-kernel-3-x86_64/gcta64 \
        --bfile /well/PROCARDIS/hrc_1kgp3/merge_p3_hrc/procardis_maf1_4gcta/procardis.maf1.chr$i \
        --chr $i \
        --cojo-file $sumstat_filepath  \
        --cojo-slct \
        --cojo-p 5e-8 \
        --cojo-collinear 0.9 \
        --cojo-wind 10000 \
        --out "$output_basename"_chr"$i" \
        --thread-num 8
done