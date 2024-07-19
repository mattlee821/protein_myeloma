#!/bin/bash

#SBATCH --job-name=extract-exposure-125k
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5-10:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

cd ~/001_projects/protein_myeloma/
  
Rscript scripts/005_colocalization/001_exposure-window/125k/001_extract-window.R