rm(list=ls())
set.seed(821)

# environment ====
library(data.table)
library(TwoSampleMR)
library(dplyr)
library(tidyr)

# data ====
data_exposure <- fread("analysis/005_colocalization/data_exposure/125k/data_exposure.txt")
data_outcome <- fread("analysis/005_colocalization/data_outcome/125k/data_outcome1.txt")

# harmonise 1====
data_harmonise <- harmonise_data(data_exposure, data_outcome, action = 2)
data_harmonise$remove_duplicates <- paste0(data_harmonise$SNP, "_", data_harmonise$id.exposure)
data_harmonise <- data_harmonise[!duplicated(data_harmonise$remove_duplicates),]
write.table(data_harmonise, "analysis/005_colocalization/data_harmonise/125k/data_harmonise-outcome1.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
