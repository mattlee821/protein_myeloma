rm(list=ls())
set.seed(821)

# environment ====
library(data.table)
library(TwoSampleMR)
library(dplyr)
library(tidyr)

# exposure data ====
LIST_FILES_UKB <- c(list.files(path = "/data/GWAS_data/work/UKB_PPP/cis-snps/combined/", pattern = "gz", all.files = T, full.names = T, recursive = T),
                    list.files(path = "/data/GWAS_data/work/UKB_PPP/cis-snps/european/", pattern = "gz", all.files = T, full.names = T, recursive = T))
LIST_FILES_DECODE <- c(list.files(path = "/data/GWAS_data/work/ferkingstad_2021_PMID34857953/cis-snps/window/", pattern = "gz", all.files = T, full.names = T, recursive = T))

## function to read and process each UKB file ====
process_UKB_file <- function(FILE) {
  # read
  data <- fread(FILE, sep = " ")
  # split column
  data <- separate(data, V16, sep = "\t", remove = TRUE, 
                   into = c("V16", "V17", "V18", "V19", "V20", "V21"))
  # column names
  colnames(data) <- c("CHR", "POS", "ID", "other_allele.exposure", "effect_allele.exposure", "eaf.exposure", 
                      "INFO", "samplesize.exposure", "TEST", "beta.exposure", "se.exposure", 
                      "CHISQ", "LOG10P", "EXTRA", "exposure", "ID", "REF", "ALT", "SNP", "POS19", "POS38")
  # remove ID col
  data <- select(data, -ID)
  # convert pval
  data$pval.exposure <- 10^(-data$LOG10P)
  # extract population from LIST_FILES
  population <- ifelse(grepl("combined", FILE), "combined", "european")
  data$population <- population
  # extract window size from LIST_FILES
  window_size <- ifelse(grepl("125k", FILE), "125k",
                        ifelse(grepl("250k", FILE), "250k",
                               ifelse(grepl("500k", FILE), "500k",
                                      ifelse(grepl("1mb", FILE), "1mb", NA))))
  data$window_size <- window_size
  # add filename
  FILE_NAME <- basename(FILE)
  FILE_NAME <- sub(".gz.", ";", FILE_NAME)
  data$id.exposure <- FILE_NAME
  # add SNP_lead column
  data$SNP_lead <- sub("^.*?;", "", data$id.exposure)
  # make id.exposure
  data$id.exposure <- paste0(data$id.exposure, ";", data$window, ";", data$population)
  # return
  return(data)
}

## function to read and process each DECODE file ====
process_DECODE_file <- function(FILE) {
  # read
  data <- fread(FILE, col.names = c("CHR", "POS", "SNPID", "SNP", "EA", "OA",
                                        "beta.exposure", "pval.exposure", "minus_log10_pval", "se.exposure", "samplesize.exposure",
                                        "EAF", "exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure"))
  # format CHR
  data$CHR <- sub("chr", "", data$CHR)
  data$CHR[data$CHR == "X"] <- "23"
  data$CHR <- as.integer(data$CHR)
  # Modify the "exposure" column
  data$exposure <- gsub("/data/protein_GWAS_ferkingstad_EU_2021/files/processed/", "", data$exposure)
  data$exposure <- gsub(".txt.gz.unzipped", "", data$exposure)
  # population
  data$population <- "DECODE"
  # extract window size from LIST_FILES
  window_size <- ifelse(grepl("125k", FILE), "125k",
                        ifelse(grepl("250k", FILE), "250k",
                               ifelse(grepl("500k", FILE), "500k",
                                      ifelse(grepl("1mb", FILE), "1mb", NA))))
  data$window_size <- window_size
  # add filename
  FILENAME <- basename(FILE)
  FILENAME <- sub(".txt.gz.annotated.gz.exclusions.gz.alleles.gz.", ";", FILENAME)
  data$id.exposure <- FILENAME
  # add SNP_lead column
  data$SNP_lead <- sub("^.*?;", "", data$id.exposure)
  # make id.exposure
  data$id.exposure <- paste0(data$id.exposure, ";", data$window, ";", data$population)
  # return
  return(data)
}

# read data ====
exposure_list1 <- lapply(LIST_FILES_UKB, process_UKB_file)
exposure_list2 <- lapply(LIST_FILES_DECODE, process_DECODE_file)
exposure_list <- c(exposure_list1, exposure_list2)
data_exposure <- bind_rows(exposure_list)
write.table(data_exposure, "analysis/005_colocalization/data_exposure.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
snps <- unique(data_exposure$SNP)

# outcome data1 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_ieu-b-4957.txt.gz"
data_outcome1 <- read_outcome_data(DATA,
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
data_outcome1$id.outcome <- data_outcome1$outcome
## transform beta from BOLT LMM scale
data_outcome1$u <- 601/372617
data_outcome1$beta.outcome <-  data_outcome1$beta.outcome/(data_outcome1$u * (1- data_outcome1$u))
data_outcome1$se.outcome <-  data_outcome1$se.outcome/(data_outcome1$u * (1- data_outcome1$u)) 

# outcome data2 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_finngen.txt.gz"
data_outcome2 <- read_outcome_data(DATA,
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
data_outcome2$id.outcome <- data_outcome2$outcome

# outcome data3 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_UKB-finngen.txt.gz"
data_outcome3 <- read_outcome_data(DATA,
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
data_outcome3$id.outcome <- data_outcome3$outcome

# outcome data4 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_UKB.txt.gz"
data_outcome4 <- read_outcome_data(DATA,
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
data_outcome4$id.outcome <- data_outcome4$outcome

# combine ====
data_outcome1$chr.outcome <- as.character(data_outcome1$chr.outcome)
data_outcome2$chr.outcome <- as.character(data_outcome2$chr.outcome)
data_outcome3$chr.outcome <- as.character(data_outcome3$chr.outcome)
data_outcome4$chr.outcome <- as.character(data_outcome4$chr.outcome)
data_outcome <- bind_rows(data_outcome1, data_outcome2, data_outcome3, data_outcome4)
write.table(data_outcome, "analysis/005_colocalization/data_outcome.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# harmonise ====
data_harmonise <- harmonise_data(data_exposure, data_outcome, action = 2)
data_harmonise$remove_duplicates <- paste0(data_harmonise$SNP, "_", data_harmonise$id.exposure)
data_harmonise <- data_harmonise[!duplicated(data_harmonise$remove_duplicates),]
write.table(data_harmonise, "analysis/005_colocalization/data_harmonise.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
