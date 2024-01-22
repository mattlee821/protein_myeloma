#!/bin/bash

#SBATCH --job-name=MR-myeloma-protein
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5-10:0:00
#SBATCH --mem=100000M

cd ~/001_projects/protein_myeloma/
Rscript scripts/004_myeloma-protein/001_MR.R          
