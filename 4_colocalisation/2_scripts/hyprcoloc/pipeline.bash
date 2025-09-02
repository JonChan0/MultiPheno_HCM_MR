#!/bin/bash

## munge input for hyprcoloc
sbatch --array=2-32 run_munge_mr.sh

## run hyprcoloc
for i in {2..32}; do bash run_hyprcoloc.bash $i; done
