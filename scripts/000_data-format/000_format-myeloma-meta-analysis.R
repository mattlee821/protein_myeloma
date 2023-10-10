rm(list=ls())
set.seed(821)

# environment ====
library(data.table)
library(dplyr)

# data ====
data <- fread("/data/GWAS_data/files/myeloma_GWAS/Multiple_myeloma/Meta_analysis_Multiple_myeloma_cancer_GWAS_formatted.csv")
data <- data[grep("^rs", data$MarkerName), ]
data <- data[grepl("\\+\\+|\\-\\-", data$Direction), ]

write.table(data, "data/myeloma_data/Meta_analysis_Multiple_myeloma_cancer_GWAS_formatted.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
