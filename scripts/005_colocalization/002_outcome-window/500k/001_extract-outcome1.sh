#!/bin/bash

#SBATCH --job-name=window-outcome1-500k
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-10:0:00
#SBATCH --mem=50000M
#SBATCH --partition=low_p

cd ~/001_projects/protein_myeloma/

Rscript scripts/005_colocalization/002_outcome-window/500k/001_extract-outcome1.R
