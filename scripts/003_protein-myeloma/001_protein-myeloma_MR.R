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
methods <- methods[c(1,3,6,10,13),1]

## exposure ====
data_exposure <- fread("analysis/001_instruments/instruments-proteins.txt")
data_exposure$f_stats <- (data_exposure$beta.exposure / data_exposure$se.exposure)^2 
write.table(data_exposure, "analysis/003_protein_myeloma/data_exposure.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# outcome
data_outcome <- fread("analysis/002_outcomes/instruments-proteins_outcome-cancer.txt")
table(data_outcome$id.outcome)
write.table(data_outcome, "analysis/003_protein_myeloma/data_outcome.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## harmonize data ====
data_harmonise <- harmonise_data(data_exposure, data_outcome, action = 2)
data_harmonise <- subset(data_harmonise, SNP_index == 1) # drop duplicate SNP from proxy search
write.table(data_harmonise, "analysis/003_protein_myeloma/data_harmonise.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## MR ====
data_mr <- mr(data_harmonise, method_list = methods)
write.table(data_mr, "analysis/003_protein_myeloma/data_mr.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## steiger ====
### data_outcome1 - myeloma-ieu-b-4957
data <- subset(data_harmonise, outcome == "myeloma-ieu-b-4957")
data$r.outcome <- get_r_from_lor(lor = data$beta.outcome, af = data$eaf.outcome, 
                                           ncase = 601, ncontrol = 372617, prevalence = (25.98/100000)*100, model = "logit")
data$r.exposure <- get_r_from_bsen(b = data$beta.exposure, se = data$se.exposure,
                                             n = min(data_exposure$samplesize.exposure))
steiger_results1 <- directionality_test(data)

### data_outcome2 - meta-analysis
data <- subset(data_harmonise, outcome == "myeloma-UKB-finngen")
data$r.outcome <- get_r_from_lor(lor = data$beta.outcome, af = data$eaf.outcome, 
                                           ncase = 1649, ncontrol = 727247, 
                                           prevalence = summary(metaprop(event = c(21.48, 25.98),
                                                                         n = c(100000, 100000), 
                                                                         studlab = c("finngen", "UKB"),
                                                                         method = "Inverse",
                                                                         sm = "PRAW"))[["TE.random"]], 
                                           model = "logit")
data$r.exposure <- get_r_from_bsen(b = data$beta.exposure, se = data$se.exposure,
                                             n = min(data_exposure$samplesize.exposure))
steiger_results2 <- directionality_test(data)

### data_outcome3 - finngen
data <- subset(data_harmonise, outcome == "myeloma-finngen")
data$r.outcome <- get_r_from_lor(lor = data$beta.outcome, af = data$eaf.outcome, 
                                           ncase = 1085, ncontrol = 271463, prevalence = (21.48/100000)*100, model = "logit")
data$r.exposure <- get_r_from_bsen(b = data$beta.exposure, se = data$se.exposure,
                                             n = min(data_exposure$samplesize.exposure))
steiger_results3 <- directionality_test(data)

### data_outcome4 - UKB
data <- subset(data_harmonise, outcome == "myeloma-UKB")
data$r.outcome <- get_r_from_lor(lor = data$beta.outcome, af = data$eaf.outcome, 
                                           ncase = 564, ncontrol = 455784, prevalence = (25.98/100000)*100, model = "logit")
data$r.exposure <- get_r_from_bsen(b = data$beta.exposure, se = data$se.exposure,
                                             n = min(data_exposure$samplesize.exposure))
steiger_results4 <- directionality_test(data)

## combine
steiger_results <- bind_rows(steiger_results1, steiger_results2, steiger_results3, steiger_results4)
write.table(steiger_results, "analysis/003_protein_myeloma/data_steiger.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
