rm(list=ls())
set.seed(821)

# environment ====
library(data.table)
library(TwoSampleMR)
library(dplyr)
library(tidyr)

# data ====
data_exposure <- fread("analysis/005_colocalization/data_exposure.txt")
snps <- fread("analysis/005_colocalization/snplist-exposure.txt")
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
data_outcome$id.outcome <- data_outcome$outcome
## transform beta from BOLT LMM scale
data_outcome$u <- 601/372617
data_outcome$beta.outcome <-  data_outcome$beta.outcome/(data_outcome$u * (1- data_outcome$u))
data_outcome$se.outcome <-  data_outcome$se.outcome/(data_outcome$u * (1- data_outcome$u)) 
write.table(data_outcome, "analysis/005_colocalization/data_outcome/data_outcome1.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
## harmonise ====
data_harmonise <- harmonise_data(data_exposure, data_outcome, action = 2)
data_harmonise$remove_duplicates <- paste0(data_harmonise$SNP, "_", data_harmonise$id.exposure)
data_harmonise <- data_harmonise[!duplicated(data_harmonise$remove_duplicates),]
write.table(data_harmonise, "analysis/005_colocalization/data_harmonise/data_harmonise1.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

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
data_outcome$id.outcome <- data_outcome$outcome
write.table(data_outcome, "analysis/005_colocalization/data_outcome/data_outcome2.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
## harmonise ====
data_harmonise <- harmonise_data(data_exposure, data_outcome, action = 2)
data_harmonise$remove_duplicates <- paste0(data_harmonise$SNP, "_", data_harmonise$id.exposure)
data_harmonise <- data_harmonise[!duplicated(data_harmonise$remove_duplicates),]
write.table(data_harmonise, "analysis/005_colocalization/data_harmonise/data_harmonise3.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# outcome data3 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_UKB-finngen.txt.gz"
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
data_outcome$id.outcome <- data_outcome$outcome
write.table(data_outcome, "analysis/005_colocalization/data_outcome/data_outcome2.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
## harmonise ====
data_harmonise <- harmonise_data(data_exposure, data_outcome, action = 2)
data_harmonise$remove_duplicates <- paste0(data_harmonise$SNP, "_", data_harmonise$id.exposure)
data_harmonise <- data_harmonise[!duplicated(data_harmonise$remove_duplicates),]
write.table(data_harmonise, "analysis/005_colocalization/data_harmonise/data_harmonise3.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# outcome data4 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_UKB.txt.gz"
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
data_outcome$id.outcome <- data_outcome$outcome
write.table(data_outcome, "analysis/005_colocalization/data_outcome/data_outcome4.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
## harmonise ====
data_harmonise <- harmonise_data(data_exposure, data_outcome, action = 2)
data_harmonise$remove_duplicates <- paste0(data_harmonise$SNP, "_", data_harmonise$id.exposure)
data_harmonise <- data_harmonise[!duplicated(data_harmonise$remove_duplicates),]
write.table(data_harmonise, "analysis/005_colocalization/data_harmonise/data_harmonise4.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

