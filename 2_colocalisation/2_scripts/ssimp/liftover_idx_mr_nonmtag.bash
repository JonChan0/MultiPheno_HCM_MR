#!/bin/bash 

idx=$1
#file=$2

chr=`grep -w ^$idx nonMTAG_indexes.txt | cut -f 2`
st=`grep -w ^$idx nonMTAG_indexes.txt | cut -f 3`
en=`grep -w ^$idx nonMTAG_indexes.txt | cut -f 4`
file=`grep -w ^$idx nonMTAG_indexes.txt | cut -f 5`

echo $chr $st $en $file

# chr=`grep -w ^$idx ../hcm_coord_b37.txt | cut -f 3`
# st=`grep -w ^$idx ../hcm_coord_b37.txt | cut -f 4`
# en=`grep -w ^$idx ../hcm_coord_b37.txt | cut -f 5`

awk -v OFS='\t' '{print "chr" $1,$2-1,$2,$5}' \
    chunks_mr_nonMTAG/${file}.chr${chr}.${st}.${en} \
    | sed '1d' \
	  > liftover_mr_nonMTAG/${file}.chr${chr}.${st}.${en}.b37.bed

/well/PROCARDIS/bin/liftover/liftOver \
    liftover_mr_nonMTAG/${file}.chr${chr}.${st}.${en}.b37.bed \
    /well/PROCARDIS/bin/liftover/hg19ToHg38.over.chain.gz \
    liftover_mr_nonMTAG/${file}.chr${chr}.${st}.${en}.b38.bed \
    liftover_mr_nonMTAG/${file}.chr${chr}.${st}.${en}.unlifted.bed

perl munge_ssimp.pl \
     chunks_mr_nonMTAG/${file}.chr${chr}.${st}.${en} \
     liftover_mr_nonMTAG/${file}.chr${chr}.${st}.${en}.b38.bed 0.3 \
     > chunks_mr_nonMTAG_output/${file}.chr${chr}.${st}.${en}.final_inf0.3.tsv \
