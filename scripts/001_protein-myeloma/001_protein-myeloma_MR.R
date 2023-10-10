rm(list=ls())
set.seed(821)
# MR analysis of measures of adiposity and endometrial cancer 

# environment ====
## library ====
#remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(data.table)
library(RadialMR)
library(dplyr)
library(tidyverse)
library(EpiViz)

### colours
#install.packages("wesanderson")
library(wesanderson)
d1 <- wes_palette("Royal1", type = "discrete")
d2 <- wes_palette("GrandBudapest2", type = "discrete")
d3 <- wes_palette("Cavalcanti1", type = "discrete")
d4 <- wes_palette("Rushmore1", type = "discrete")
discrete_wes_pal <- c(d1, d2, d3, d4)
rm(d1,d2,d3,d4)

## extract exposure instruments ====
data <- read_exposure_data("data/protein_data/cis_snps.txt",
                           clump = F,
                           sep = "\t",
                           snp_col = "rsID",
                           beta_col = "BETA",
                           se_col = "SE",
                           eaf_col = "EAF",
                           effect_allele_col = "EffectAllele",
                           other_allele_col = "OtherAllele",
                           pval_col = "Pval",
                           samplesize_col = "N",
                           phenotype = "phenotype",
                           min_pval = 5e-9)
data <- data[!duplicated(data[,"exposure"]),]

# sequence_id column 
ncols <- max(stringr::str_count(data$exposure, "_")) + 1
colmn <- paste0("col", 1:ncols)
data <-
  tidyr::separate(
    data = data,
    col = exposure,
    sep = "_",
    into = colmn,
    remove = FALSE)
data$sequence_id <- paste0(data$col1, "_", data$col2)
# remove extra cols 
data <- data[, -which(names(data) %in% c("col1","col2","col3","col4","col5","col6","col7","col8","col9","col10"))]

# exposure data ====
data$exposure <- as.factor(data$exposure)
data$f_stats <- (data$beta.exposure / data$se.exposure)^2 
exposure_data <- data
write.table(exposure_data, "analysis/001_protein_myeloma/exposure_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# outcome data ====
## extract outcome data ====
outcome_data1 <- read_outcome_data("data/myeloma_data/ieu-b-4957.tsv",
  snps = exposure_data$SNP,
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  min_pval = 1e-200,
  log_pval = FALSE,
  chr_col = "chromosome",
  pos_col = "base_pair_location")
outcome_data1$outcome <- "ieu-b-4957"
outcome_data1$id.outcome <- "ieu-b-4957"
## transform beta from BLOT LMM scale
outcome_data1$u <- 601/372617
outcome_data1$beta.outcome <-  outcome_data1$beta.outcome/(outcome_data1$u * (1- outcome_data1$u))
outcome_data1$se.outcome <-  outcome_data1$se.outcome/(outcome_data1$u * (1- outcome_data1$u)) 

outcome_data2 <- read_outcome_data("data/myeloma_data/Meta_analysis_Multiple_myeloma_cancer_GWAS_formatted.txt",
                                  snps = exposure_data$SNP,
                                  sep = "\t",
                                  snp_col = "MarkerName",
                                  beta_col = "beta",
                                  se_col = "se",
                                  eaf_col = "Freq1",
                                  effect_allele_col = "Allele1",
                                  other_allele_col = "Allele2",
                                  pval_col = "P-value",
                                  min_pval = 1e-200,
                                  log_pval = FALSE)
outcome_data2$outcome <- "myeloma_meta-analysis"
outcome_data2$id.outcome <- "myeloma_meta-analysis"

outcome_data <- bind_rows(outcome_data1, outcome_data2)
write.table(outcome_data, "analysis/001_protein_myeloma/outcome_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## harmonize data ====
harmonise_data <- harmonise_data(exposure_data, outcome_data, action = 2)
write.table(harmonise_data, "analysis/001_protein_myeloma/harmonise_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## MR ====
mr_results <- mr(harmonise_data)
write.table(mr_results, "analysis/001_protein_myeloma/mr_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## steiger ====
### outcome_data1
harmonise_data <- harmonise_data(exposure_data, outcome_data1, action = 2)
harmonise_data$r.outcome <- get_r_from_lor(lor = harmonise_data$beta.outcome, af = harmonise_data$eaf.outcome, 
                                           ncase = 601, ncontrol = 372617, prevalence = (24000/70000000)*100, model = "logit")
harmonise_data$r.exposure <- get_r_from_bsen(b = harmonise_data$beta.exposure, se = harmonise_data$se.exposure,
                                             n = min(exposure_data$samplesize.exposure))
steiger_results <- directionality_test(harmonise_data)
write.table(steiger_results, "analysis/001_protein_myeloma/steiger_results_ieu-b-4957.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

### outcome_data2
harmonise_data <- harmonise_data(exposure_data, outcome_data2, action = 2)
harmonise_data$r.outcome <- get_r_from_lor(lor = harmonise_data$beta.outcome, af = harmonise_data$eaf.outcome, 
                                           ncase = 1649, ncontrol = 727247, prevalence = (24000/70000000)*100, model = "logit")
harmonise_data$r.exposure <- get_r_from_bsen(b = harmonise_data$beta.exposure, se = harmonise_data$se.exposure,
                                             n = min(exposure_data$samplesize.exposure))
steiger_results <- directionality_test(harmonise_data)
write.table(steiger_results, "analysis/001_protein_myeloma/steiger_results_myeloma-meta-analysis.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# format ====
data <- mr_results
## OR and CI
data$OR <- exp(data$b)
data$lower_ci <- exp(data$b - (1.96 * data$se))
data$upper_ci <- exp(data$b + (1.96 * data$se))

## sequence_id column 
ncols <- max(stringr::str_count(data$exposure, "_")) + 1
colmn <- paste0("col", 1:ncols)
data <-
  tidyr::separate(
    data = data,
    col = exposure,
    sep = "_",
    into = colmn,
    remove = FALSE)
data$sequence_id <- paste0(data$col1, "_", data$col2)

## protein_name column 
data$id.exposure <- sub(".*?_", "", data$exposure) # remove first part of sequence_id from column
data$id.exposure <- sub(".*?_", "", data$id.exposure) # remove second part of sequence_id from column
data$protein <- sub(".*?_", "", data$id.exposure) # remove genename

## gene name column 
data$gene <- gsub("\\_.*","",data$id.exposure)

## remove extra cols 
data <- select(data, !colmn)

## organise data frame
data <- data[,c("exposure", "outcome", "sequence_id", "protein", "gene", "nsnp", "method", "b", "se", "pval", "OR", "lower_ci", "upper_ci")]

## save
write.table(data, "analysis/001_protein_myeloma/001_MR_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
