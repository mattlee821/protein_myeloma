rm(list=ls())

# environment ====
library(genetics.binaRies)
library(plinkbinr)
library(coloc)
library(viridis)
library(data.table)
library(ieugwasr)
library(dplyr)
library(TwoSampleMR)
library(tidyverse)
source("scripts/003_colocalization/002_coloc/functions/my_coloc_chriswallace.R")

# data ====
files <- list.files("analysis/003_colocalization/results", pattern = "harmonise_data", recursive = TRUE, full.names = TRUE)
harmonise_data_list <- fread("analysis/003_colocalization/results/harmonise_data.txt")
harmonise_data_list <- split(harmonise_data_list, harmonise_data_list$id.exposure)
tempdir()

# coloc ====
table_master <- data.frame() # make empty dataframe for final results

for (i in 1245:length(harmonise_data_list)){
  
  label_exposure <- unique(harmonise_data_list[[i]]$exposure)
  label <- paste0(label_exposure, "_", "myeloma")
  label_outcome <- "myeloma"  
  # make ld matrix ====
  ld <- ld_matrix_local(
    harmonise_data_list[[i]]$SNP,
    with_alleles = FALSE, 
    bfile = "/data/GWAS_data/files/references/1kG_v3/EUR/EUR",
    plink_bin = get_plink_exe())
  # Identify rows and columns with no missing values
  non_missing_rows <- apply(ld, 1, function(x) !all(is.nan(x)))
  non_missing_cols <- apply(ld, 2, function(x) !all(is.nan(x)))
  ld <- ld[non_missing_rows, non_missing_cols]
  
  # format LD matrix and harmonised list ====
  ld <- ld[which(rownames(ld) %in% harmonise_data_list[[i]]$SNP), which(colnames(ld) %in% harmonise_data_list[[i]]$SNP)]
  harmonise_data_list[[i]] <- harmonise_data_list[[i]][which(harmonise_data_list[[i]]$SNP %in% rownames(ld)),]
  ld <- ld[match(harmonise_data_list[[i]]$SNP,rownames(ld)),]
  ld <- ld[,match(harmonise_data_list[[i]]$SNP, colnames(ld))]
  harmonise_data_list[[i]] <- harmonise_data_list[[i]][match(rownames(ld), harmonise_data_list[[i]]$SNP),]
  
  # make lists for coloc ====
  coloc_data_exposure <- list(beta = harmonise_data_list[[i]]$beta.exposure, varbeta = harmonise_data_list[[i]]$se.exposure^2, MAF = harmonise_data_list[[i]]$eaf.exposure, type = "quant", N = 35559, snp = rownames(ld), LD = ld, position = harmonise_data_list[[i]]$POS)
  coloc_data_outcome <- list(beta = harmonise_data_list[[i]]$beta.outcome, varbeta = harmonise_data_list[[i]]$se.outcome^2, MAF = harmonise_data_list[[i]]$eaf.outcome, type = "cc", N = 120328, snp = rownames(ld), LD = ld, position = harmonise_data_list[[i]]$POS)
  
  # coloc ====  
  coloc_results <- coloc.abf(dataset1 = coloc_data_exposure, dataset2 = coloc_data_outcome, p1 = 1E-6, p2 = 1E-6, p12 = 1E-7)
  
  pdf(paste0("analysis/003_colocalization/results/figures/", label, ".pdf"), 
      height = 10, width = 10)
  coloc_sensitivity <- my_sensitivity(coloc_results, "H4 > 0.9", 
                                      trait1_title = label_exposure, trait2_title = label_outcome)
  dev.off()
  
  # save ====
  saveRDS(coloc_results, paste0("analysis/003_colocalization/results/", label, ".RData"))
  
  # make table ====
  table <- data.frame(
    exposure = label_exposure,
    outcome = label_outcome,
    id = label,
    nsnps = coloc_results["summary"][[1]][1],
    h0 = coloc_results["summary"][[1]][2],
    h1 = coloc_results["summary"][[1]][3],
    h2 = coloc_results["summary"][[1]][4],
    h3 = coloc_results["summary"][[1]][5],
    h4 = coloc_results["summary"][[1]][6])
  
  table_master <- rbind(table_master, table)
  
}

write.table(table_master, "analysis/003_colocalization/results/001_coloc_results_missing.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
