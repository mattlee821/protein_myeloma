#!/bin/bash

#SBATCH --job-name=filelist-45
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

# set directory ====
export SNP_LIST=~/001_projects/protein_myeloma/analysis/001_instruments/snplist-cancer.txt
export OUTCOMES=~/001_projects/protein_myeloma/analysis/002_outcomes/filelist/

cd ${OUTCOMES}

tmp=$(mktemp) || { ret="$?"; printf 'Failed to create temp file\n'; exit "$ret"; }
VAR1=filelist-45
  
for file in $(cat "$VAR1"); do
  echo "running $(basename "$file")"
  zgrep -w -F -f "$SNP_LIST" "$file" > "$tmp" &&
  mv -- "$tmp" "${OUTCOMES}$(basename "$file")"
  echo "finished $(basename "$file")"
done
