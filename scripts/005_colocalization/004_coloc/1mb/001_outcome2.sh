#!/bin/bash

#SBATCH --job-name=coloc-1mb-2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10-10:0:00
#SBATCH --mem=60000M
#SBATCH --partition=low_p

export TMPDIR=/scratch/leem/temp/ # this will change temp location
cd /data/MET_share/work/001_projects/protein_myeloma/
Rscript scripts/005_colocalization/004_coloc/1mb/001_outcome2.R              
