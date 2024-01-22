#!/bin/bash

#SBATCH --job-name=extract_cis_snps
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M

cd /data/protein_GWAS_ferkingstad_EU_2021/files/

mkdir cis_snps_1mb/

export SNP_FILE=/data/protein_GWAS_ferkingstad_EU_2021/work/cis_snp_list.txt

export FILE_LIST=/data/protein_GWAS_ferkingstad_EU_2021/work/filelist_cis_snps.txt

for file in `cat $FILE_LIST`; do
    
    gzip -d -c ${file} > ${file}.unzipped

    awk -v diff=1000000 '

function abs(x) { return (x < 0.0) ? -x : x }

BEGIN         { FS=OFS=" " }

FNR==NR       { if (FNR>1)                      # 1st file; skip header and ...
                   bp_list[$7][$1]=$2           # save contents in our bp_list[FILE][CHR] array
                next
              }

FNR==1        { close(outfile)                  # close previous output file
                fn=FILENAME
                outfile=fn ".cis"
                if (fn in bp_list)              # if fn in 1st file then ...
                   print > outfile              # print header else ...
                else                            # skip to next input file; also addresses gwas* matching on gwas*_extract files, ie, these will be skipped, too
                   nextfile
                next
              }

fn in bp_list { if ($1 in bp_list[fn] && abs(bp_list[fn][$1] - $2) <= diff)
                   print > outfile
              }
' $SNP_FILE ${file}.unzipped

    mv ${file}.unzipped.cis cis_snps_1mb/

    rm ${file}.unzipped 
done

cd /data/protein_GWAS_ferkingstad_EU_2021/files/cis_snps_1mb

ls *unzipped.cis > filelist

for file in `cat filelist`; do
    tail -n +2 ${file} > ${file}.txt
    rm ${file}
done