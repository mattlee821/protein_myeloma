#!/bin/bash

#SBATCH --job-name=make-multiple-filelists
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

# set directory ====
cd ~/001_projects/protein_myeloma/analysis/002_outcomes/
mkdir filelist
cd filelist

# make file list combining ukb and ferkingstad protein GWAS files ====
rm filelist*
ls /data/GWAS_data/work/UKB_PPP/european/*gz > filelist-ukb-eu
ls /data/GWAS_data/work/UKB_PPP/combined/*gz > filelist-ukb-combined
ls /data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS/*gz > filelist-ferkingstad
wc -l filelist-ukb-eu
wc -l filelist-ukb-combined
wc -l filelist-ferkingstad
cat filelist-ukb-eu filelist-ukb-combined filelist-ferkingstad > filelist
wc -l filelist
rm filelist-ukb-eu filelist-ukb-combined filelist-ferkingstad

# set number of protein files to split ====
# adjust "lines_per_file" based on number of protein files, "wc -l filelist-ferkingstad" + "wc -l filelist-ukb"
lines_per_file=$((($(wc -l < filelist) / 99) + 1))
echo $lines_per_file

# split to create the new files  ====
input_file="filelist"
output_prefix="filelist-"
split --lines="$lines_per_file" --numeric-suffixes=1 "$input_file" "$output_prefix"
wc -l *
rm filelist
