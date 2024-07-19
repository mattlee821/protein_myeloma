rm(list=ls())

outcome <- "ieu-b-4957"
path <- "analysis/003_colocalization/results/ieu-b-4957/"
filenames <- list.files(path = paste0(path, "/"), pattern = ".RData", full.names = TRUE)

exposure <- gsub(".RData", "", filenames)
exposure <- gsub("analysis/003_colocalization/results/ieu-b-4957///", "", exposure)
exposure <- gsub(paste0("_",outcome), "", exposure)

table_master <- data.frame() # make empty dataframe for final results

for (i in 1:length(filenames)){
  coloc_results <- readRDS(filenames[i])
  
  label_exposure <- exposure[i]
  label_outcome <- outcome
  label <- paste0(label_exposure, "_", label_outcome)
  
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

table_outcome1 <- table_master

outcome <- "myeloma-meta-analysis"
path <- "analysis/003_colocalization/results/myeloma-meta-analysis/"
filenames <- list.files(path = paste0(path, "/"), pattern = ".RData", full.names = TRUE)

exposure <- gsub(".RData", "", filenames)
exposure <- gsub("analysis/003_colocalization/results/myeloma-meta-analysis///", "", exposure)
exposure <- gsub(paste0("_",outcome), "", exposure)

table_master <- data.frame() # make empty dataframe for final results

for (i in 1:length(filenames)){
  coloc_results <- readRDS(filenames[i])
  
  label_exposure <- exposure[i]
  label_outcome <- outcome
  label <- paste0(label_exposure, "_", label_outcome)
  
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

table_outcome2 <- table_master

# combined 
table <- rbind(table_outcome1, table_outcome2)

write.table(table, "analysis/003_colocalization/results/001_coloc_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

