rm(list=ls())
set.seed(821)

# environment ====
#remotes::install_github("MRCIEU/TwoSampleMR")
# remotes::install_github("mattlee821/functions", force = T)
library(functions)
library(TwoSampleMR)
library(data.table)
library(RadialMR)
library(dplyr)
library(tidyverse)
library(EpiViz)
library(meta)
library(wesanderson)

### methods
methods <- mr_method_list()
methods_heterogeneity <- subset(methods, heterogeneity_test == TRUE)$obj
methods_heterogeneity <- methods_heterogeneity[c(1,2,3,5)]
methods <- methods[c(3,6,10,13),1]

# exposure data ====
data_exposure <- read.table("analysis/001_instruments/instruments-cancer.txt", header = T, sep = "\t")
write.table(data_exposure, "analysis/004_myeloma_protein/data_exposure.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# outcome data ====
## decode proteins
data <- fread("analysis/002_outcomes/instruments-cancer_outcome-proteinsDECODE.txt", col.names = c(
  "chr.outcome",
  "pos.outcome",
  "name",
  "SNP",
  "EA_incorrect",
  "OA_incorrect",
  "beta.outcome",
  "pval.outcome",
  "min log10P",
  "se.outcome",
  "samplesize.outcome",
  "impMAF",
  "outcome",
  "effect_allele.outcome",
  "other_allele.outcome",
  "eaf.outcome"
))
data$chr.outcome <- gsub(pattern = "chr", "", data$chr.outcome)
data$chr.outcome <- as.numeric(data$chr.outcome)
data1 <- select(data, chr.outcome, pos.outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, SNP, beta.outcome, se.outcome, pval.outcome, samplesize.outcome, outcome)
data1$outcome <- paste0(data1$outcome, ";DECODE")

## UKB proteins eu
data <- fread("analysis/002_outcomes/instruments-cancer_outcome-proteinsUKB-eu.txt", sep = " ", col.names = c(
  "chr.outcome",
  "pos.outcome",
  "ID",
  "other_allele.outcome",
  "effect_allele.outcome",
  "eaf.outcome",
  "INFO",
  "samplesize.outcome",
  "TEST",
  "beta.outcome",
  "se.outcome",
  "CHISQ",
  "LOG10P",
  "EXTRA",
  "phenotype",
  "V1"
))
data <- separate(data, V1, into = c("SNPID", "REF", "ALT", "SNP", "POS19", "POS38"), sep = "\t")
data <- TwoSampleMR::format_data(data, type = "outcome", snps = NULL,  header = TRUE,  
                    phenotype_col = "phenotype",  
                    snp_col = "SNP",  
                    beta_col = "beta.outcome",  
                    se_col = "se.outcome",  
                    eaf_col = "eaf.outcome",  
                    effect_allele_col = "effect_allele.outcome",  
                    other_allele_col = "other_allele.outcome",  
                    pval_col = "LOG10P", 
                    id_col = "phenotype", 
                    chr_col = "chr.outcome",  
                    pos_col = "pos.outcome",  
                    samplesize_col = "samplesize.outcome",
                    log_pval = T)
data2 <- select(data, chr.outcome, pos.outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, SNP, beta.outcome, se.outcome, pval.outcome, samplesize.outcome, outcome)
data2$outcome <- paste0(data2$outcome, ";UKB-EU")

## UKB proteins combined
data <- fread("analysis/002_outcomes/instruments-cancer_outcome-proteinsUKB-combined.txt", sep = " ", col.names = c(
  "chr.outcome",
  "pos.outcome",
  "ID",
  "other_allele.outcome",
  "effect_allele.outcome",
  "eaf.outcome",
  "INFO",
  "samplesize.outcome",
  "TEST",
  "beta.outcome",
  "se.outcome",
  "CHISQ",
  "LOG10P",
  "EXTRA",
  "phenotype",
  "V1"
))
data <- separate(data, V1, into = c("SNPID", "REF", "ALT", "SNP", "POS19", "POS38"), sep = "\t")
data <- TwoSampleMR::format_data(data, type = "outcome", snps = NULL,  header = TRUE,  
                    phenotype_col = "phenotype",  
                    snp_col = "SNP",  
                    beta_col = "beta.outcome",  
                    se_col = "se.outcome",  
                    eaf_col = "eaf.outcome",  
                    effect_allele_col = "effect_allele.outcome",  
                    other_allele_col = "other_allele.outcome",  
                    pval_col = "LOG10P", 
                    id_col = "phenotype", 
                    chr_col = "chr.outcome",  
                    pos_col = "pos.outcome",  
                    samplesize_col = "samplesize.outcome",
                    log_pval = T)
data3 <- select(data, chr.outcome, pos.outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, SNP, beta.outcome, se.outcome, pval.outcome, samplesize.outcome, outcome)
data3$outcome <- paste0(data3$outcome, ";UKB-EU")

## combine and format
data_outcome <- bind_rows(data1, data2, data3)
data_outcome$outcome <- gsub("/data/protein_GWAS_ferkingstad_EU_2021/files/processed/", "", data_outcome$outcome)
data_outcome$outcome <- gsub(".txt.gz.unzipped", "", data_outcome$outcome)
data_outcome$id.outcome <- data_outcome$outcome
data_outcome1 <- data_outcome %>% # remap the EU/ALL label based on samplesize
  group_by(SNP) %>%
  mutate(
    outcome = case_when(
      samplesize.outcome == max(samplesize.outcome) ~ sub("UKB-ALL|UKB-EU", "UKB-ALL", outcome),
      samplesize.outcome == min(samplesize.outcome) ~ sub("UKB-ALL|UKB-EU", "UKB-EU", outcome),
      TRUE ~ outcome
    ),
    id.outcome = outcome
  ) %>%
  ungroup()
write.table(data_outcome, "analysis/004_myeloma_protein/data_outcome.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# # harmonize data ====
# data_harmonise <- harmonise_data(data_exposure, data_outcome, action = 2)
# data_harmonise <- data_harmonise %>%
#   filter(mr_keep == TRUE)
# write.table(data_harmonise, "analysis/004_myeloma_protein/data_harmonise.txt", 
#             row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
# 
# MR ====
data_harmonise <- fread("analysis/004_myeloma_protein/data_harmonise.txt")
data_mr <- mr(data_harmonise)
write.table(data_mr, "analysis/004_myeloma_protein/data_mr.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# steiger ====
## data_outcome1 - ieu
data <- subset(data_harmonise, exposure == "myeloma-ieu-b-4957")
data$r.exposure <- get_r_from_lor(lor = data$beta.exposure, af = data$eaf.exposure, 
                                   ncase = 601, ncontrol = 372617, prevalence = (25.98/100000)*100, model = "logit")
data$r.outcome <- get_r_from_bsen(b = data$beta.outcome, se = data$se.outcome,
                                   n = min(data$samplesize.outcome))
steiger_results1 <- directionality_test(data)

## data_outcome2 - meta-analysis
data <- subset(data_harmonise, exposure == "myeloma-UKB-finngen")
data$r.exposure <- get_r_from_lor(lor = data$beta.exposure, af = data$eaf.exposure, 
                                   ncase = 1649, ncontrol = 727247, 
                                   prevalence = summary(metaprop(event = c(21.48, 25.98),
                                                                 n = c(100000, 100000), 
                                                                 studlab = c("finngen", "UKB"),
                                                                 method = "Inverse",
                                                                 sm = "PRAW"))[["TE.random"]], 
                                   model = "logit")
data$r.outcome <- get_r_from_bsen(b = data$beta.outcome, se = data$se.outcome,
                                   n = min(data$samplesize.outcome))
steiger_results2 <- directionality_test(data)

## data_outcome3 - finngen
data <- subset(data_harmonise, exposure == "myeloma-finngen")
data$r.exposure <- get_r_from_lor(lor = data$beta.exposure, af = data$eaf.exposure, 
                                   ncase = 1085, ncontrol = 271463, prevalence = (21.48/100000)*100, model = "logit")
data$r.outcome <- get_r_from_bsen(b = data$beta.outcome, se = data$se.outcome,
                                   n = min(data_harmonise$samplesize.outcome))
steiger_results3 <- directionality_test(data)

## data_outcome4 - UKB
data <- subset(data_harmonise, exposure == "myeloma-UKB")
data$r.exposure <- get_r_from_lor(lor = data$beta.exposure, af = data$eaf.exposure, 
                                   ncase = 564, ncontrol = 455784, prevalence = (25.98/100000)*100, model = "logit")
data$r.outcome <- get_r_from_bsen(b = data$beta.outcome, se = data$se.outcome,
                                   n = min(data_harmonise$samplesize.outcome))
steiger_results4 <- directionality_test(data)
## combine
steiger_results <- bind_rows(steiger_results1, steiger_results2, steiger_results3, steiger_results4)
write.table(steiger_results, "analysis/004_myeloma_protein/data_steiger.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
