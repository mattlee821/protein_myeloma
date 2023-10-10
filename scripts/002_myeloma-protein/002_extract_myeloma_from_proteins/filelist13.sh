#!/bin/bash

#SBATCH --job-name=filelist13
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M

# extract
export SNP_LIST=~/001_projects/protein_myeloma/analysis/002_myeloma_protein/snplist.txt 
export OUTCOMES=~/001_projects/protein_myeloma/analysis/002_myeloma_protein/outcome_data/

cd /data/protein_GWAS_ferkingstad_EU_2021/files/processed/

tmp=$(mktemp) || { ret="$?"; printf 'Failed to create temp file\n'; exit "$ret"; }
for file in `cat filelist13`; do
	zgrep -w -F -f ${SNP_LIST} ${file} > ${tmp} &&
	mv -- ${tmp} ${file}
    mv ${file} ${OUTCOMES}
done
