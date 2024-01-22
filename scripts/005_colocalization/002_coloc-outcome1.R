rm(list=ls())
set.seed(821)

# environment ====
# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# remotes::install_github("MRCIEU/genetics.binaRies", force = F)
# remotes::install_github("explodecomputer/plinkbinr", force = F)
# remotes::install_github("chr1swallace/coloc@main", force = F)
# remotes::install_github("sjmgarnier/viridis", force = F)
library(genetics.binaRies)
library(plinkbinr)
library(coloc)
library(viridis)
library(data.table)
library(ieugwasr)
library(dplyr)
library(TwoSampleMR)
library(tidyverse)

source("scripts/005_colocalization/002_coloc/functions/my_coloc_chriswallace.R")

# data ====
data_harmonise <- fread("analysis/005_colocalization/data_harmonise.txt")
outcome <- "myeloma-ieu-b-4957"
data_harmonise <- subset(data_harmonise, id.outcome == outcome)
list_data_harmonise <- split(data_harmonise, data_harmonise$id.exposure)

# loop over all harmonised data and run ld matrix, formatting, coloc, save ====
table_master <- data.frame() # make empty dataframe for final results

for (i in 1:length(list_data_harmonise)) {

  # variables ====
  label <- paste0(unique(list_data_harmonise[[i]]$exposure), ";", unique(list_data_harmonise[[i]]$SNP_lead), ";", outcome)
  label_exposure <- paste0(unique(list_data_harmonise[[i]]$exposure), ";", unique(list_data_harmonise[[i]]$SNP_lead))
  window_size <- unique(list_data_harmonise[[i]]$window_size)
  population <- unique(list_data_harmonise[[i]]$population)
  
  # make a dataframe for empty data (sometimes there are no variants for coloc and so this will give us an NA row for those situations) ====
  row_info <- data.frame(
    exposure = label_exposure,
    outcome = outcome,
    id = label,
    nsnps = NA,
    h0 = NA,
    h1 = NA,
    h2 = NA,
    h3 = NA,
    h4 = NA,
    window_size = window_size,
    population = population
  )
  
  # LD matrix (do inside tryCatch() as sometime syou have no variants) ====
  tryCatch({
    # make ld matrix ====
    ld <- ld_matrix_local(
      list_data_harmonise[[i]]$SNP,
      with_alleles = FALSE, 
      bfile = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR",
      plink_bin = get_plink_exe())
    
    # check if no variants remain after --extract
    if (dim(ld)[1] == 0) {
      # Handle case where no variants remain
      cat("Warning: No variants remaining after --extract for file", i, "\n")
    } else {
      
      # format LD matrix and harmonised list ====
      ld <- ld[which(rownames(ld) %in% list_data_harmonise[[i]]$SNP), which(colnames(ld) %in% list_data_harmonise[[i]]$SNP)]
      list_data_harmonise[[i]] <- list_data_harmonise[[i]][which(list_data_harmonise[[i]]$SNP %in% rownames(ld)),]
      ld <- ld[match(list_data_harmonise[[i]]$SNP,rownames(ld)),]
      ld <- ld[,match(list_data_harmonise[[i]]$SNP, colnames(ld))]
      list_data_harmonise[[i]] <- list_data_harmonise[[i]][match(rownames(ld), list_data_harmonise[[i]]$SNP),]
      
      # make lists for coloc ====
      coloc_data_exposure <- list(beta = list_data_harmonise[[i]]$beta.exposure, varbeta = list_data_harmonise[[i]]$se.exposure^2, MAF = list_data_harmonise[[i]]$eaf.exposure, type = "quant", N = 35559, snp = rownames(ld), LD = ld, position = list_data_harmonise[[i]]$POS)
      coloc_data_outcome <- list(beta = list_data_harmonise[[i]]$beta.outcome, varbeta = list_data_harmonise[[i]]$se.outcome^2, MAF = list_data_harmonise[[i]]$eaf.outcome, type = "cc", N = 120328, snp = rownames(ld), LD = ld, position = list_data_harmonise[[i]]$POS)
      
      # coloc ====  
      coloc_results <- coloc.abf(dataset1 = coloc_data_exposure, dataset2 = coloc_data_outcome, p1 = 1E-6, p2 = 1E-6, p12 = 1E-7)
      # coloc sensitivity 
      pdf(paste0("analysis/005_colocalization/results/test/", population, "/", window_size, "/figures/", label, ".pdf"), 
          height = 10, width = 10)
      coloc_sensitivity <- my_sensitivity(coloc_results, "H4 > 0.8", 
                                          trait1_title = label_exposure, trait2_title = outcome)
      dev.off()
      # save
      saveRDS(coloc_results, paste0("analysis/005_colocalization/results/test/", population, "/", window_size, "/", label, ".RData"))
      
      # make results table =====
      row_info$nsnps <- coloc_results["summary"][[1]][1]
      row_info$h0 <- coloc_results["summary"][[1]][2]
      row_info$h1 <- coloc_results["summary"][[1]][3]
      row_info$h2 <- coloc_results["summary"][[1]][4]
      row_info$h3 <- coloc_results["summary"][[1]][5]
      row_info$h4 <- coloc_results["summary"][[1]][6]
    }
    
  }, error = function(e) {
    cat("Error occurred in file", i, ":", conditionMessage(e), "\n")
  })
  
  # add results table (row_info) to table_master (outside of the tryCatch block)
  table_master <- rbind(table_master, row_info)
}

# write
write.table(table_master, paste0("analysis/005_colocalization/results/test/001_coloc-results", outcome, ".txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


