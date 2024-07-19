#!/bin/bash

#SBATCH --job-name=harm-outcome4-125k
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5-10:0:00
#SBATCH --mem=50000M
#SBATCH --partition=low_p

cd ~/001_projects/protein_myeloma/

Rscript scripts/005_colocalization/003_harmonise-window/125k/001_harmonise-outcome4.R
