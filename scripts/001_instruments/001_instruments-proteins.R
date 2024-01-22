rm(list=ls())
set.seed(821)

# environment ====
library(TwoSampleMR)
library(dplyr)

# data ====
LIST_FILES <- c(list.files(path = "/data/GWAS_data/work/UKB_PPP/cis-snps/", pattern = "002", all.files = T, full.names = T),
                list.files(path = "/data/GWAS_data/work/ferkingstad_2021_PMID34857953/cis-snps/", pattern = "004_cissnps-exposure.txt", all.files = T, full.names = T))
LIST_FILES <- LIST_FILES[!grepl("nonEU", LIST_FILES)]
list_data <- list()

# file list ====
for (FILE in LIST_FILES){
  
  # make label
  label <- gsub("/data/GWAS_data/work/UKB_PPP/cis-snps//002_cissnps-", "", FILE)
  label <- gsub("/data/GWAS_data/work/ferkingstad_2021_PMID34857953/cis-snps//", "", label)
  label <- gsub(".txt", "", label)
  # Map label
  label <- case_when(
    label == "combined-allancestries" ~ "UKB-allancestries",
    label == "discovery-EU" ~ "UKB-EU",
    label == "replicationn-nonEU" ~ "UKB-nonEU",
    label == "004_cissnps-exposure" ~ "DECODE",
    TRUE ~ label  # Default to original label if not matched
  )
  
  # data
  data <- read_exposure_data(filename = FILE, sep = "\t", 
                             clump = F, 
                             min_pval = 5e-08,
                             phenotype_col = "phenotype", 
                             id_col = "phenotype", 
                             snp_col = "SNP", 
                             beta_col = "BETA", 
                             se_col = "SE", 
                             pval_col = "P", log_pval = TRUE,
                             eaf_col = "EAF", 
                             effect_allele_col = "EA", 
                             other_allele_col = "OA", 
                             chr_col = "CHR", 
                             pos_col = "POS",
                             samplesize_col = "N")
  data$id.exposure <- paste0(data$exposure, ";", label)
  list_data <- append(list_data, list(data))
}

# combine ====
data <- bind_rows(list_data)
data <- data[!is.na(data$other_allele.exposure), ]
write.table(data, "analysis/001_instruments/instruments-proteins.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data <- data %>% 
  distinct(SNP)
write.table(data, "analysis/001_instruments/snplist-proteins.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
