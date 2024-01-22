#!/bin/bash

#SBATCH --job-name=extract-outcome
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-10:0:00
#SBATCH --mem=500000M
#SBATCH --partition=low_p

cd ~/001_projects/protein_myeloma/

Rscript scripts/005_colocalization/002_extract-outcome.R
