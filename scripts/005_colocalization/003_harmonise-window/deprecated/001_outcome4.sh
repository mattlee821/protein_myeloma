#!/bin/bash

#SBATCH --job-name=harmonise-outcome4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5-10:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

cd ~/001_projects/protein_myeloma/

Rscript scripts/005_colocalization/003_harmonise-window/001_outcome4.R
