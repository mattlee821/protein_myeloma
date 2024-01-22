rm(list=ls())
set.seed(821)

# environment ====
#remotes::install_github("MRCIEU/TwoSampleMR")
# remotes::install_github("mattlee821/functions", force = T)
library(functions)
library(TwoSampleMR)
library(data.table)
library(ieugwasr)
library(gwasvcf)
library(genetics.binaRies)
gwasvcf::set_plink(genetics.binaRies::get_plink_binary())

## SNPs ====
data_exposure <- fread("analysis/001_instruments/snplist-proteins.txt", header = F, col.names = "SNP")

# outcome data1 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_ieu-b-4957.txt.gz"
data_outcome <- read_outcome_data(DATA,
                                   snps = data_exposure$SNP,
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
## extract proxies
data_outcome <- functions::proxy_search(data_exposure = data_exposure, 
                                        data_outcome = data_outcome, 
                                        data_outcome_path = DATA, 
                                        data_reference = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR.bim", 
                                        data_reference_path = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR", 
                                        tag_r2 = 0.8, 
                                        tag_kb = 5000, 
                                        tag_nsnp = 5000)
## transform beta from BOLT LMM scale
data_outcome$u <- 601/372617
data_outcome$beta.outcome <-  data_outcome$beta.outcome/(data_outcome$u * (1- data_outcome$u))
data_outcome$se.outcome <-  data_outcome$se.outcome/(data_outcome$u * (1- data_outcome$u)) 
##
data_outcome1 <- data_outcome
data_outcome1$id.outcome <- data_outcome1$outcome

# outcome data2 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_finngen.txt.gz"
data_outcome <- read_outcome_data(DATA,
                                  snps = data_exposure$SNP,
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
## extract proxies
data_outcome2 <- functions::proxy_search(data_exposure = data_exposure, 
                                        data_outcome = data_outcome, 
                                        data_outcome_path = DATA, 
                                        data_reference = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR.bim", 
                                        data_reference_path = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR", 
                                        tag_r2 = 0.8, 
                                        tag_kb = 5000, 
                                        tag_nsnp = 5000)
data_outcome2$id.outcome <- data_outcome2$outcome

# outcome data3 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_UKB-finngen.txt.gz"
data_outcome <- read_outcome_data(DATA,
                                  snps = data_exposure$SNP,
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
## extract proxies
data_outcome3 <- functions::proxy_search(data_exposure = data_exposure, 
                                         data_outcome = data_outcome, 
                                         data_outcome_path = DATA, 
                                         data_reference = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR.bim", 
                                         data_reference_path = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR", 
                                         tag_r2 = 0.8, 
                                         tag_kb = 5000, 
                                         tag_nsnp = 5000)
data_outcome3$id.outcome <- data_outcome3$outcome

# outcome data4 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_UKB.txt.gz"
data_outcome <- read_outcome_data(DATA,
                                  snps = data_exposure$SNP,
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
## extract proxies
data_outcome4 <- functions::proxy_search(data_exposure = data_exposure, 
                                         data_outcome = data_outcome, 
                                         data_outcome_path = DATA, 
                                         data_reference = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR.bim", 
                                         data_reference_path = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR", 
                                         tag_r2 = 0.8, 
                                         tag_kb = 5000, 
                                         tag_nsnp = 5000)
data_outcome4$id.outcome <- data_outcome4$outcome

# combine ====
data_outcome1$chr.outcome <- as.character(data_outcome1$chr.outcome)
data_outcome2$chr.outcome <- as.character(data_outcome2$chr.outcome)
data_outcome3$chr.outcome <- as.character(data_outcome3$chr.outcome)
data_outcome4$chr.outcome <- as.character(data_outcome4$chr.outcome)
data_outcome <- bind_rows(data_outcome1, data_outcome2, data_outcome3, data_outcome4)
write.table(data_outcome, "analysis/002_outcomes/instruments-proteins_outcome-cancer.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
