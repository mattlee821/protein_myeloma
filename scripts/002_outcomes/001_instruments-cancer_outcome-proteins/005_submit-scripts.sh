#!/bin/bash

#SBATCH --job-name=make-multiple-submissionscripts
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

# set directory ====
cd ~/001_projects/protein_myeloma/scripts/002_outcomes/001_instruments-cancer_outcome-proteins/004_filelist

# submit the filelist-*.sh scripts  ====
for file in filelist-*.sh; do
sbatch $file
sleep 2
done