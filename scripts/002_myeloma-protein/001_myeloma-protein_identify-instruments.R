rm(list=ls())
set.seed(821)

# environment ====
## library ====
#remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(data.table)
library(RadialMR)
library(dplyr)

## extract exposure instruments ====
## exposure 1
exposure_data1 <- read_exposure_data("data/myeloma_data/ieu-b-4957.tsv",
  clump = F,
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  chr_col = "chromosome",
  pos_col = "base_pair_location")
exposure_data1$exposure <- "ieu-b-4957"
exposure_data1$id.exposure <- "ieu-b-4957"

### subset
a <- subset(exposure_data1, pval.exposure <= 5e-8)
b <- subset(exposure_data1, pval.exposure <= 5e-7)
c <- subset(exposure_data1, pval.exposure <= 5e-6)
nrow(a) # 0
nrow(b) # 3
nrow(c) # 63

###
exposure_data <- subset(exposure_data1, pval.exposure <= 5e-7)
exposure_data$f_stats <- (exposure_data$beta.exposure / exposure_data$se.exposure)^2 
exposure_data %>%
  group_by(exposure) %>%
  summarise(mean = mean(f_stats))
exposure_data_clumped <- clump_data(exposure_data, clump_kb = 10000, clump_r2 = 0.001) # removed 1 SNP for absence from reference (not because of LD)
exposure_data1 <- exposure_data

## exposure 2
exposure_data2 <- read_exposure_data("data/myeloma_data/Meta_analysis_Multiple_myeloma_cancer_GWAS_formatted.txt",
                                     clump = F,
                                     sep = "\t",
                                     snp_col = "MarkerName",
                                     beta_col = "beta",
                                     se_col = "se",
                                     eaf_col = "Freq1",
                                     effect_allele_col = "Allele1",
                                     other_allele_col = "Allele2",
                                     pval_col = "P-value",
                                     chr_col = "chromosome")
exposure_data2$exposure <- "myeloma_meta-analysis"
exposure_data2$id.exposure <- "myeloma_meta-analysis"

### subset
a <- subset(exposure_data2, pval.exposure <= 5e-8)
b <- subset(exposure_data2, pval.exposure <= 5e-7)
c <- subset(exposure_data2, pval.exposure <= 5e-6)
nrow(a) # 5
nrow(b) # 26
nrow(c) # 133

###
exposure_data <- subset(exposure_data2, pval.exposure <= 5e-8)
exposure_data$f_stats <- (exposure_data$beta.exposure / exposure_data$se.exposure)^2 
exposure_data %>%
  group_by(exposure) %>%
  summarise(mean = mean(f_stats))
exposure_data_clumped <- clump_data(exposure_data, clump_kb = 10000, clump_r2 = 0.001) # removed 4 SNPs for LD

# manually add chr and pos from dbSNP
exposure_data_clumped$chr.exposure <- 7
exposure_data_clumped$pos.exposure <- 21939032
exposure_data2 <- exposure_data_clumped

## combine
exposure_data <- rbind(exposure_data1, exposure_data2)
write.table(exposure_data, "analysis/002_myeloma_protein/exposure_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## make snplist 
snps <- as.data.frame(exposure_data$SNP)
write.table(snps, "analysis/002_myeloma_protein/snplist.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
