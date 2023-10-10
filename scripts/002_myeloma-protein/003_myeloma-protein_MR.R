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
exposure_data <- read.table("analysis/002_myeloma_protein/exposure_data.txt", header = T, sep = "\t")
### transform form bolt lmm scale
exposure_data$u <- 601/372617
exposure_data$beta.exposure[exposure_data$exposure == "ieu-b-4957"] <- exposure_data$beta.exposure[exposure_data$exposure == "ieu-b-4957"] / (exposure_data$u * (1 - exposure_data$u))[exposure_data$exposure == "ieu-b-4957"]
exposure_data$se.exposure[exposure_data$exposure == "ieu-b-4957"] <- exposure_data$se.exposure[exposure_data$exposure == "ieu-b-4957"] / (exposure_data$u * (1 - exposure_data$u))[exposure_data$exposure == "ieu-b-4957"]

# outcome data ====
## extract outcome data ====
data <- read.table("analysis/002_myeloma_protein/outcome_data/outcome_data.txt", col.names = c(
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
data <- select(data, chr.outcome, pos.outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, SNP, beta.outcome, se.outcome, pval.outcome, samplesize.outcome, outcome)

## sequence_id column 
data$outcome <- gsub("/data/protein_GWAS_ferkingstad_EU_2021/files/processed/", "", data$outcome)
data$outcome <- gsub(".txt.gz.unzipped", "", data$outcome)
ncols <- max(stringr::str_count(data$outcome, "_")) + 1
colmn <- paste0("col", 1:ncols)
data <-
  tidyr::separate(
    data = data,
    col = outcome,
    sep = "_",
    into = colmn,
    remove = FALSE)
data$sequence_id <- paste0(data$col1, "_", data$col2)
# remove extra cols 
data <- data[, -which(names(data) %in% c("col1","col2","col3","col4","col5","col6","col7","col8","col9","col10"))]
data$id.outcome <- data$outcome
outcome_data <- data
write.table(outcome_data, "analysis/002_myeloma_protein/outcome_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## harmonize data ====
harmonise_data <- harmonise_data(exposure_data, outcome_data, action = 2)
write.table(harmonise_data, "analysis/002_myeloma_protein/harmonise_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## MR ====
mr_results <- mr(harmonise_data)
write.table(mr_results, "analysis/002_myeloma_protein/mr_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## steiger ====
### exposure_data1
exposure_data1 <- subset(exposure_data, exposure == "ieu-b-4957")
harmonise_data <- harmonise_data(exposure_data1, outcome_data, action = 2)
harmonise_data$r.exposure <- get_r_from_lor(lor = harmonise_data$beta.outcome, af = harmonise_data$eaf.outcome, 
                                           ncase = 601, ncontrol = 372617, prevalence = (24000/70000000)*100, model = "logit")
harmonise_data$r.outcome <- get_r_from_bsen(b = harmonise_data$beta.outcome, se = harmonise_data$se.outcome,
                                             n = min(outcome_data$samplesize.outcome))
steiger_results <- directionality_test(harmonise_data)
write.table(steiger_results, "analysis/002_myeloma_protein/steiger_results_ieu-b-4957.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

### exposure_data2
exposure_data1 <- subset(exposure_data, exposure == "myeloma_meta-analysis")
harmonise_data <- harmonise_data(exposure_data1, outcome_data, action = 2)
harmonise_data$r.exposure <- get_r_from_lor(lor = harmonise_data$beta.outcome, af = harmonise_data$eaf.outcome, 
                                            ncase = 1649, ncontrol = 727247, prevalence = (24000/70000000)*100, model = "logit")
harmonise_data$r.outcome <- get_r_from_bsen(b = harmonise_data$beta.outcome, se = harmonise_data$se.outcome,
                                            n = min(outcome_data$samplesize.outcome))
steiger_results <- directionality_test(harmonise_data)
write.table(steiger_results, "analysis/002_myeloma_protein/steiger_results_myeloma-meta-analysis.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# format ====
data <- mr_results
## OR and CI
data$lower_ci <- data$b - (1.96 * data$se)
data$upper_ci <- data$b + (1.96 * data$se)

## sequence_id column 
ncols <- max(stringr::str_count(data$outcome, "_")) + 1
colmn <- paste0("col", 1:ncols)
data <-
  tidyr::separate(
    data = data,
    col = outcome,
    sep = "_",
    into = colmn,
    remove = FALSE)
data$sequence_id <- paste0(data$col1, "_", data$col2)

## protein_name column 
data$id.outcome <- sub(".*?_", "", data$outcome) # remove first part of sequence_id from column
data$id.outcome <- sub(".*?_", "", data$id.outcome) # remove second part of sequence_id from column
data$protein <- sub(".*?_", "", data$id.outcome) # remove genename

## gene name column 
data$gene <- gsub("\\_.*","",data$id.outcome)

## remove extra cols 
data <- select(data, !colmn)

## organise data frame
data <- data[,c("exposure", "outcome", "sequence_id", "protein", "gene", "nsnp", "method", "b", "se", "pval", "lower_ci", "upper_ci")]

## save
write.table(data, "analysis/002_myeloma_protein/001_MR_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
