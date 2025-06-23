#!/bin/bash

idx=$1;
chr=`grep -w ^${idx} hcm_nonMTAG_b38.txt1 | cut -f 3`;
st=`grep -w ^${idx} hcm_nonMTAG_b38.txt1 | cut -f 4`;
en=`grep -w ^${idx} hcm_nonMTAG_b38.txt1 | cut -f 5`;
loci=`grep -w ^${idx} hcm_nonMTAG_b38.txt1 | cut -f 6`;

echo $idx $chr $st $en $loci
perl ../scripts/munge_hyprcoloc_input_h.pl $chr $st $en datasets_nonMTAG.txt ${loci}MR_nonMTAG

