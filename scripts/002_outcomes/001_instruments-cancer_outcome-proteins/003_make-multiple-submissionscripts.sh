#!/bin/bash

#SBATCH --job-name=make-multiple-submissionscripts
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

# set directory ====
cd ~/001_projects/protein_myeloma/scripts/002_outcomes/001_instruments-cancer_outcome-proteins
mkdir 004_filelist
cd 004_filelist
rm filelist*
rm filenames
  
# create a file with the names of all files within a directory  ====
ls ~/001_projects/protein_myeloma/analysis/002_outcomes/filelist/ > filenames

# create multiple .sh scripts from a single file with names in based on a master script  ====
# line (NR == 17) holds the variable name
# line (NR == 3) holds the job name
cat filenames | while read i; do echo ${i}; 
awk '{ 
  if (NR == 3) 
    print "#SBATCH --job-name='${i}'";
  else if (NR == 17)
    print "VAR1='${i}'";
  else
    print $0
}' ../002_master-filelist.sh > ${i}.sh; done 
rm filenames