#!/bin/bash

#SBATCH --job-name=coloc-125k-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10-10:0:00
#SBATCH --mem=60000M
#SBATCH --partition=low_p

export TMPDIR=/scratch/leem/temp/ # this will change temp location
cd /data/MET_share/work/001_projects/protein_myeloma/
Rscript scripts/005_colocalization/004_coloc/125k/001_outcome1.R              
