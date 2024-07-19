#!/bin/bash

#SBATCH --job-name=harmonise-outcome3new
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5-10:0:00
#SBATCH --mem=150000M
#SBATCH --partition=low_p

cd ~/001_projects/protein_myeloma/

Rscript scripts/005_colocalization/003_harmonise-window/001_outcome3.R
