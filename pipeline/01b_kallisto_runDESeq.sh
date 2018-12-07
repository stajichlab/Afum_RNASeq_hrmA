#!/bin/bash
#SBATCH --nodes 1 --ntasks 1 --out logs/DEseq.log --mem 4G  --time 2:00:00 -p short

Rscript Rscripts/kallisto_profile.R
