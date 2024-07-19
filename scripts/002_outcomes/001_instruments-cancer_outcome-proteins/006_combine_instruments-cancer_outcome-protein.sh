FILE_IN=~/001_projects/protein_myeloma/analysis/002_outcomes/filelist/
FILE_OUT=~/001_projects/protein_myeloma/analysis/002_outcomes/

cat ${FILE_IN}*european* > ${FILE_OUT}instruments-cancer_outcome-proteinsUKB-eu.txt
cat ${FILE_IN}*combined* > ${FILE_OUT}instruments-cancer_outcome-proteinsUKB-combined.txt
cat ${FILE_IN}*DECODE* > ${FILE_OUT}instruments-cancer_outcome-proteinsDECODE.txt
