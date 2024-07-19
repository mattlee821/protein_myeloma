#!/bin/bash

#SBATCH --job-name=filelist-99
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

# set directory ====
export SNP_LIST=~/001_projects/protein_myeloma/analysis/001_instruments/snplist-cancer.txt
export OUTCOMES=~/001_projects/protein_myeloma/analysis/002_outcomes/filelist/
cd "${OUTCOMES}" || exit

tmp=$(mktemp) || { ret="$?"; printf 'Failed to create temp file\n'; exit "$ret"; }
VAR1=filelist-99

for filelist in $VAR1; do
  while IFS= read -r file; do
    echo "running $(basename "$file")"
    directory=$(dirname "$file")
    if [[ $directory == *"european"* ]]; then
      output_name="${OUTCOMES}$(basename "$file")_european.txt"
    elif [[ $directory == *"combined"* ]]; then
      output_name="${OUTCOMES}$(basename "$file")_combined.txt"
    else
      output_name="${OUTCOMES}$(basename "$file")_DECODE.txt"
    fi
    zgrep -w -F -f "$SNP_LIST" "$file" > "$tmp" &&
    mv -- "$tmp" "$output_name"
    echo "finished $(basename "$file")"
  done < "$filelist"
done
