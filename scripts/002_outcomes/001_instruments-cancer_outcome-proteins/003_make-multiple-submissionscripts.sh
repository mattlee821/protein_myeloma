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
rm *

# create a file with the names of all files within a directory  ====
ls ~/001_projects/protein_myeloma/analysis/002_outcomes/filelist/filelist* > filenames

# create multiple .sh scripts from a single file with names in based on a master script  ====
# line (NR == 17) holds the variable name
# line (NR == 3) holds the job name
cat filenames | while read i; do
    echo ${i}; 
    filename=$(basename ${i})
    awk '{
        if (NR == 3)
            print "#SBATCH --job-name='"${filename}"'";
        else if (NR == 16)
            print "VAR1='"${filename}"'";
        else
            print $0
    }' ../002_master-filelist.sh > ~/001_projects/protein_myeloma/scripts/002_outcomes/001_instruments-cancer_outcome-proteins/004_filelist/${filename}.sh
done 
rm filenames
