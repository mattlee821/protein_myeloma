rm(list=ls())
set.seed(821)

library(dplyr)
library(data.table)

LIST_FILES <- c(
  list.files(path = "/data/MET_share/work/001_projects/protein_myeloma/analysis/005_colocalization/data_results/125k/", pattern = "001_coloc", all.files = T, full.names = T, recursive = T),
  list.files(path = "/data/MET_share/work/001_projects/protein_myeloma/analysis/005_colocalization/data_results/250k/", pattern = "001_coloc", all.files = T, full.names = T, recursive = T),
  list.files(path = "/data/MET_share/work/001_projects/protein_myeloma/analysis/005_colocalization/data_results/500k/", pattern = "001_coloc", all.files = T, full.names = T, recursive = T),
  list.files(path = "/data/MET_share/work/001_projects/protein_myeloma/analysis/005_colocalization/data_results/1mb/", pattern = "001_coloc", all.files = T, full.names = T, recursive = T)
)

list_data <- lapply(LIST_FILES, fread, header = T, sep = "\t")
data <- bind_rows(list_data)
write.table(data, "/data/MET_share/work/001_projects/protein_myeloma/analysis/005_colocalization/data_results/001_results-coloc.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

