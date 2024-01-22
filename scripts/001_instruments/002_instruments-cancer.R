rm(list=ls())
set.seed(821)

# environment ====
## library ====
#remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(data.table)
library(RadialMR)
library(dplyr)
library(ieugwasr)

## extract exposure instruments ====
## exposure 1 ====
data_exposure1 <- read_exposure_data("/data/GWAS_data/work/myeloma/myeloma_combined_ieu-b-4957.txt.gz",
                                     clump = F,
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
                                     phenotype_col = "phenotype", id_col = "phenotype")
data_exposure1$exposure <- data_exposure1$exposure

### subset
a <- subset(data_exposure1, pval.exposure <= 5e-8)
b <- subset(data_exposure1, pval.exposure <= 5e-7)
c <- subset(data_exposure1, pval.exposure <= 5e-6)
nrow(a) # 0
nrow(b) # 3
nrow(c) # 63
###
data_exposure <- subset(data_exposure1, pval.exposure <= 5e-7)
data_exposure$f_stats <- (data_exposure$beta.exposure / data_exposure$se.exposure)^2 
data_exposure %>%
  group_by(exposure) %>%
  summarise(mean = mean(f_stats))
### clump
data_exposure_clumped <- data_exposure
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "SNP"] <- "rsid"
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "pval.exposure"] <- "pval"
data_exposure_clumped <- ld_clump(dat = data_exposure_clumped,
                         clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8,
                         pop = "EUR",
                         access_token = NULL,
                         bfile = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR",
                         plink_bin = genetics.binaRies::get_plink_binary())
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "rsid"] <- "SNP"
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "pval"] <- "pval.exposure"
data_exposure1 <- data_exposure # use unclumped data

## exposure 2 ====
data_exposure2 <- read_exposure_data("/data/GWAS_data/work/myeloma/myeloma_combined_UKB-finngen.txt.gz",
                                     clump = F,
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
                                     phenotype_col = "phenotype", id_col = "phenotype")
data_exposure2$exposure <- data_exposure2$exposure

### subset
a <- subset(data_exposure2, pval.exposure <= 5e-8)
b <- subset(data_exposure2, pval.exposure <= 5e-7)
c <- subset(data_exposure2, pval.exposure <= 5e-6)
nrow(a) # 5
nrow(b) # 26
nrow(c) # 133

###
data_exposure <- subset(data_exposure2, pval.exposure <= 5e-8)
data_exposure$f_stats <- (data_exposure$beta.exposure / data_exposure$se.exposure)^2 
data_exposure %>%
  group_by(exposure) %>%
  summarise(mean = mean(f_stats))
### clump
data_exposure_clumped <- data_exposure
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "SNP"] <- "rsid"
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "pval.exposure"] <- "pval"
data_exposure_clumped <- ld_clump(dat = data_exposure_clumped,
                                  clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8,
                                  pop = "EUR",
                                  access_token = NULL,
                                  bfile = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR",
                                  plink_bin = genetics.binaRies::get_plink_binary())
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "rsid"] <- "SNP"
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "pval"] <- "pval.exposure"
data_exposure2 <- data_exposure_clumped

## exposure 3 ====
data_exposure3 <- read_exposure_data("/data/GWAS_data/work/myeloma/myeloma_combined_finngen.txt.gz",
                                     clump = F,
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
                                     phenotype_col = "phenotype", id_col = "phenotype")
data_exposure3$exposure <- data_exposure3$exposure

### subset
a <- subset(data_exposure3, pval.exposure <= 5e-8)
b <- subset(data_exposure3, pval.exposure <= 5e-7)
c <- subset(data_exposure3, pval.exposure <= 5e-6)
nrow(a) # 0
nrow(b) # 10
nrow(c) # 131
###
data_exposure <- subset(data_exposure3, pval.exposure <= 5e-7)
data_exposure$f_stats <- (data_exposure$beta.exposure / data_exposure$se.exposure)^2 
data_exposure %>%
  group_by(exposure) %>%
  summarise(mean = mean(f_stats))
### clump
data_exposure_clumped <- data_exposure
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "SNP"] <- "rsid"
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "pval.exposure"] <- "pval"
data_exposure_clumped <- ld_clump(dat = data_exposure_clumped,
                                  clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8,
                                  pop = "EUR",
                                  access_token = NULL,
                                  bfile = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR",
                                  plink_bin = genetics.binaRies::get_plink_binary())
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "rsid"] <- "SNP"
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "pval"] <- "pval.exposure"
data_exposure3 <- data_exposure %>% # this returns the single clumped SNP and the other SNPs with lowest pval that are either the only SNP on the chr or have the lowest pval on the CHR
  group_by(chr.exposure) %>%
  slice(which.min(pval.exposure)) %>%
  ungroup()

## exposure 4 ====
data_exposure4 <- read_exposure_data("/data/GWAS_data/work/myeloma/myeloma_combined_UKB.txt.gz",
                                     clump = F,
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
                                     phenotype_col = "phenotype", id_col = "phenotype")
data_exposure4$exposure <- data_exposure4$exposure

### subset
a <- subset(data_exposure4, pval.exposure <= 5e-8)
b <- subset(data_exposure4, pval.exposure <= 5e-7)
c <- subset(data_exposure4, pval.exposure <= 5e-6)
nrow(a) # 5
nrow(b) # 3
nrow(c) # 63

###
data_exposure <- subset(data_exposure4, pval.exposure <= 5e-8)
data_exposure$f_stats <- (data_exposure$beta.exposure / data_exposure$se.exposure)^2 
data_exposure %>%
  group_by(exposure) %>%
  summarise(mean = mean(f_stats))
### clump
data_exposure_clumped <- data_exposure
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "SNP"] <- "rsid"
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "pval.exposure"] <- "pval"
data_exposure_clumped <- ld_clump(dat = data_exposure_clumped,
                                  clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8,
                                  pop = "EUR",
                                  access_token = NULL,
                                  bfile = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR",
                                  plink_bin = genetics.binaRies::get_plink_binary())
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "rsid"] <- "SNP"
colnames(data_exposure_clumped)[colnames(data_exposure_clumped) == "pval"] <- "pval.exposure"
data_exposure4 <- data_exposure_clumped

## rearrange ====
data_exposure1 <- data_exposure1 %>%
  select(exposure, id.exposure, chr.exposure, pos.exposure, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure,
         beta.exposure, se.exposure, pval.exposure, f_stats)
data_exposure2 <- data_exposure2 %>%
  select(exposure, id.exposure, chr.exposure, pos.exposure, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure,
         beta.exposure, se.exposure, pval.exposure, f_stats)
data_exposure3 <- data_exposure3 %>%
  select(exposure, id.exposure, chr.exposure, pos.exposure, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure,
         beta.exposure, se.exposure, pval.exposure, f_stats)
data_exposure4 <- data_exposure4 %>%
  select(exposure, id.exposure, chr.exposure, pos.exposure, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure,
         beta.exposure, se.exposure, pval.exposure, f_stats)

# Assuming data_exposure1, data_exposure2, etc. are your data frames
data_exposure1$chr.exposure <- as.character(data_exposure1$chr.exposure)
data_exposure2$chr.exposure <- as.character(data_exposure2$chr.exposure)
data_exposure3$chr.exposure <- as.character(data_exposure3$chr.exposure)
data_exposure4$chr.exposure <- as.character(data_exposure4$chr.exposure)
# Bind rows after converting the column types
data_exposure <- bind_rows(data_exposure1, data_exposure2, data_exposure3, data_exposure4)
write.table(data_exposure, "analysis/001_instruments/instruments-cancer.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## make snplist 
data <- select(data_exposure, SNP)
write.table(data, "analysis/001_instruments/snplist-cancer.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
