rm(list=ls())
set.seed(821)

# environment ====
library(data.table)
library(TwoSampleMR)
library(dplyr)
library(tidyr)

# data ====
snps <- fread("analysis/005_colocalization/data_exposure/1mb/snplist-exposure.txt")
snps <- unique(snps$x)

# outcome data2 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_finngen.txt.gz"
data_outcome <- read_outcome_data(DATA,
                                  snps = snps,
                                  sep = "\t",
                                  snp_col = "SNP",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  eaf_col = "EAF",
                                  effect_allele_col = "EA",
                                  other_allele_col = "OA",
                                  pval_col = "P",
                                  chr_col = "CHR",
                                  pos_col = "POS", 
                                  samplesize_col = "N",  
                                  phenotype_col = "phenotype", id_col = "phenotype")
write.table(data_outcome, "analysis/005_colocalization/data_outcome/1mb/data_outcome2.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
