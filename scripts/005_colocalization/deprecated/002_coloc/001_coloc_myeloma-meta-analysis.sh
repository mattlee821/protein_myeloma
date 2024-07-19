#!/bin/bash

#SBATCH --job-name=coloc-mma
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5-10:0:00
#SBATCH --mem=100000M

cd ~/001_projects/protein_myeloma/
Rscript scripts/003_colocalization/002_coloc/001_coloc_myeloma-meta-analysis.R                
