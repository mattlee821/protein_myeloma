rm(list=ls())
set.seed(821)

# environment ====
#remotes::install_github("MRCIEU/TwoSampleMR")
# remotes::install_github("mattlee821/functions", force = T)
library(functions)
library(TwoSampleMR)
library(data.table)
library(ieugwasr)
library(gwasvcf)
library(genetics.binaRies)
gwasvcf::set_plink(genetics.binaRies::get_plink_binary())

proxy_search_test <- function(data_exposure, data_outcome, data_outcome_path, data_reference, data_reference_path,
                              tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000) {
  
  # Parameter Validation
  if (!file.exists(data_outcome_path) || !file.exists(paste0(data_reference))) {
    stop("Invalid file paths provided.")
  }
  
  if (!is.data.frame(data_exposure) || !is.data.frame(data_outcome)) {
    stop("data_exposure and data_outcome must be data frames.")
  }
  
  # exposure snps missing from outcome ====
  snps_missing <- setdiff(unique(as.factor(data_exposure$SNP)), unique(as.factor(data_outcome$SNP)))
  if (length(snps_missing) == 0) {
    message("No missing SNPs. Returning the original dataframe.")
    return(data_outcome)
  }
  
  # look-up missing SNPs in reference panel ====
  message(paste0("# 1. looking up ", length(unique(as.factor(snps_missing))), " missing-SNP(s) in the reference panel"))
  reference_header <- readLines(data_reference, n = 1) ## read the first row to identify the column containing "rs*"
  column_index <- grep("rs", strsplit(reference_header, "\t")[[1]]) ## get column index
  ## bash command to filter rows based on the identified column
  cmd <- sprintf("awk '$%d ~ /^(%s)$/' %s", column_index, paste(snps_missing, collapse = "|"), paste0(data_reference))
  ## use system to run the bash command and then read the result with fread
  reference <- data.table::fread(cmd = cmd, quote = "")
  snps_reference <- intersect(unique(as.factor(snps_missing)), unique(as.factor(reference$V2)))
  message(paste0("## ", length(unique(as.factor(snps_reference))), " of ", length(unique(as.factor(snps_missing))), " missing-SNP(s) available in the reference panel"))
  
  # find proxies for available SNPs ====
  message(paste0("# 2. extracting proxy-SNP(s) for the ", length(unique(as.factor(snps_reference))), " missing-SNP(s) from the reference panel"))
  gwasvcf::set_plink(genetics.binaRies::get_plink_binary())
  proxies <- functions::get_ld_proxies(rsid = snps_missing,
                                       bfile = data_reference_path,
                                       searchspace = NULL,
                                       tag_kb = tag_kb,
                                       tag_nsnp = tag_nsnp,
                                       tag_r2 = tag_r2,
                                       threads = 1,
                                       out = tempfile())
  ## format proxy data: change column order and names, add proxy.outcome = TRUE
  proxies <- proxies %>%
    dplyr::select(target_snp.outcome = SNP_A,
                  proxy_snp.outcome = SNP_B,
                  target_a1.outcome = A1,
                  target_a2.outcome = A2,
                  proxy_a1.outcome = B1,
                  proxy_a2.outcome = B2,
                  R) %>%
    dplyr::mutate(proxy.outcome = TRUE,
                  SNP = proxy_snp.outcome) %>%
    dplyr::select(proxy.outcome, everything())
  message(paste0("## proxy-SNP(s) for ", length(unique(proxies$target_snp.outcome)), " missing-SNP(s) found; ", "proxy-SNP(s) for ", length(unique(as.factor(snps_reference))) - length(unique(as.factor(proxies$target_snp.outcome))), " missing-SNP(s) not available (e.g., no proxy-SNP or r2 < provided)"))
  
  # extract proxies from outcome ====
  message(paste0("# 3. extracting proxy-SNP(s) from outcome"))
  proxy_snps <- proxies %>% # select unique proxy SNPs to extract
    dplyr::distinct(proxy_snp.outcome) %>%
    dplyr::pull(proxy_snp.outcome)
  data_outcome_proxies <- TwoSampleMR::read_outcome_data(filename = data_outcome_path,
                                                         snps = proxy_snps,
                                                         sep = "\t",
                                                         phenotype_col = "phenotype",
                                                         snp_col = "SNP",
                                                         beta_col = "BETA",
                                                         se_col = "SE",
                                                         eaf_col = "EAF",
                                                         effect_allele_col = "EA",
                                                         other_allele_col = "OA",
                                                         pval_col = "P",
                                                         samplesize_col = "N",
                                                         id_col = "phenotype",
                                                         chr_col = "CHR",
                                                         pos_col = "POS")
  data_outcome_proxies <- dplyr::left_join(data_outcome_proxies, proxies, by = c("SNP" = "SNP"))
  message(paste0("## proxy-SNP(s) for ", length(unique(as.factor(data_outcome_proxies$target_snp.outcome))), " of ", length(unique(as.factor(proxies$target_snp.outcome))), " missing-SNP(s) extracted"))
  
  # select proxy-SNP(s) with the highest R2 ====
  data_outcome_proxies <- data_outcome_proxies %>%
    dplyr::group_by(target_snp.outcome) %>%
    dplyr::filter(R == max(R)) %>%
    dplyr::slice(1) %>%
    dplyr::select(-R)
  
  ## Bind rows of data_outcome with data_outcome_proxies
  data_outcome <- dplyr::bind_rows(data_outcome, data_outcome_proxies)
  
  return(data_outcome)
}

## SNPs ====
data_exposure <- fread("analysis/001_instruments/instruments-proteins.txt", header = T)

# outcome data1 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_ieu-b-4957.txt.gz"
data_outcome <- read_outcome_data(DATA,
                                   snps = data_exposure$SNP,
                                   sep = "\t",
                                   snp_col = "SNP",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   eaf_col = "EAF",
                                   effect_allele_col = "EA",
                                   other_allele_col = "OA",
                                   pval_col = "P",
                                   chr_col = "CHR",
                                   pos_col = "POS", 
                                   samplesize_col = "N",  
                                   phenotype_col = "phenotype", id_col = "phenotype")
## extract proxies
data_outcome <- proxy_search_test(data_exposure = data_exposure,
                                  data_outcome = data_outcome, 
                                  data_outcome_path = DATA, 
                                  data_reference = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR.bim", 
                                  data_reference_path = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR", 
                                  tag_r2 = 0.8, 
                                  tag_kb = 5000, 
                                  tag_nsnp = 5000)
## transform beta from BOLT LMM scale
data_outcome$u <- 601/372617
data_outcome$beta.outcome <-  data_outcome$beta.outcome/(data_outcome$u * (1- data_outcome$u))
data_outcome$se.outcome <-  data_outcome$se.outcome/(data_outcome$u * (1- data_outcome$u)) 
data_outcome1 <- data_outcome

# outcome data2 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_finngen.txt.gz"
data_outcome <- read_outcome_data(DATA,
                                  snps = data_exposure$SNP,
                                  sep = "\t",
                                  snp_col = "SNP",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  eaf_col = "EAF",
                                  effect_allele_col = "EA",
                                  other_allele_col = "OA",
                                  pval_col = "P",
                                  chr_col = "CHR",
                                  pos_col = "POS", 
                                  samplesize_col = "N",  
                                  phenotype_col = "phenotype", id_col = "phenotype")
## extract proxies
data_outcome2 <- proxy_search_test(data_exposure = data_exposure, 
                                   data_outcome = data_outcome, 
                                   data_outcome_path = DATA, 
                                   data_reference = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR.bim", 
                                   data_reference_path = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR", 
                                   tag_r2 = 0.8, 
                                   tag_kb = 5000, 
                                   tag_nsnp = 5000)

# outcome data3 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_UKB-finngen.txt.gz"
data_outcome <- read_outcome_data(DATA,
                                  snps = data_exposure$SNP,
                                  sep = "\t",
                                  snp_col = "SNP",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  eaf_col = "EAF",
                                  effect_allele_col = "EA",
                                  other_allele_col = "OA",
                                  pval_col = "P",
                                  chr_col = "CHR",
                                  pos_col = "POS", 
                                  samplesize_col = "N",  
                                  phenotype_col = "phenotype", id_col = "phenotype")
## extract proxies
data_outcome3 <- proxy_search_test(data_exposure = data_exposure, 
                                   data_outcome = data_outcome, 
                                   data_outcome_path = DATA, 
                                   data_reference = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR.bim", 
                                   data_reference_path = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR", 
                                   tag_r2 = 0.8, 
                                   tag_kb = 5000, 
                                   tag_nsnp = 5000)

# outcome data4 ====
DATA <- "/data/GWAS_data/work/myeloma/myeloma_combined_UKB.txt.gz"
data_outcome <- read_outcome_data(DATA,
                                  snps = data_exposure$SNP,
                                  sep = "\t",
                                  snp_col = "SNP",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  eaf_col = "EAF",
                                  effect_allele_col = "EA",
                                  other_allele_col = "OA",
                                  pval_col = "P",
                                  chr_col = "CHR",
                                  pos_col = "POS", 
                                  samplesize_col = "N",  
                                  phenotype_col = "phenotype", id_col = "phenotype")
## extract proxies
data_outcome4 <- proxy_search_test(data_exposure = data_exposure, 
                                   data_outcome = data_outcome, 
                                   data_outcome_path = DATA, 
                                   data_reference = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR.bim", 
                                   data_reference_path = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR",
                                   tag_r2 = 0.8, 
                                   tag_kb = 5000, 
                                   tag_nsnp = 5000)

# combine ====
data_outcome1$chr.outcome <- as.character(data_outcome1$chr.outcome)
data_outcome2$chr.outcome <- as.character(data_outcome2$chr.outcome)
data_outcome3$chr.outcome <- as.character(data_outcome3$chr.outcome)
data_outcome4$chr.outcome <- as.character(data_outcome4$chr.outcome)
data_outcome <- dplyr::bind_rows(data_outcome1, data_outcome2, data_outcome3, data_outcome4)
write.table(data_outcome, "analysis/002_outcomes/instruments-proteins_outcome-cancer.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
