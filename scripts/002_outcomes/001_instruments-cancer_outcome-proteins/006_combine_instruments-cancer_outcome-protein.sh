FILE_IN=~/001_projects/protein_myeloma/analysis/002_outcomes/filelist/
FILE_OUT=~/001_projects/protein_myeloma/analysis/002_outcomes/

cat ${FILE_IN}*OID* > ${FILE_OUT}instruments-cancer_outcome-proteinsUKB.txt
cat ${FILE_IN}*.txt.gz.annotated.gz.exclusions.gz.alleles.gz > ${FILE_OUT}instruments-cancer_outcome-proteinsDECODE.txt
