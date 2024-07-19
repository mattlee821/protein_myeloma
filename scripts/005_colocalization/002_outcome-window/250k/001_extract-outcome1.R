rm(list=ls())
set.seed(821)

# environment ====
library(data.table)
library(TwoSampleMR)
library(dplyr)
library(tidyr)

# data ====
snps <- fread("analysis/005_colocalization/data_exposure/250k/snplist-exposure.txt")
snps <- unique(snps$x)

# outcome data1 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_ieu-b-4957.txt.gz"
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
## transform beta from BOLT LMM scale
data_outcome$u <- 601/372617
data_outcome$beta.outcome <-  data_outcome$beta.outcome/(data_outcome$u * (1- data_outcome$u))
data_outcome$se.outcome <-  data_outcome$se.outcome/(data_outcome$u * (1- data_outcome$u)) 
write.table(data_outcome, "analysis/005_colocalization/data_outcome/250k/data_outcome1.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")