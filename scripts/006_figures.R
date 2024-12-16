rm(list=ls())
set.seed(821)

library(ggplot2)
library(tidyverse)
library(stringr)
library(openxlsx)
library(scales)
library(dplyr)


## read in MR results
mr_results_full <- read.delim("003_protein_myeloma/data_mr.txt")

## read in harmonised data
harmonised_data <- read.delim("003_protein_myeloma/data_harmonise.txt")
harmonised_data <- harmonised_data[,c("SNP", "exposure", "outcome")]

## Restrict protein (deCODE) -> myeloma (UKB-finngen)
mr_soma_protein_myeloma_all <- subset(mr_results_full, grepl("DECODE", exposure)) 
mr_soma_protein_ukb_myeloma <- mr_soma_protein_myeloma_all[mr_soma_protein_myeloma_all$outcome == "myeloma-UKB-finngen" ,]

## Restrict protein (UKB) -> myeloma (finngen)
mr_olink_protein_myeloma_all <- mr_results_full[str_detect(mr_results_full$exposure, "^[A-Z]" ), ]
mr_olink_protein_finngen_myeloma <-  mr_olink_protein_myeloma_all[mr_olink_protein_myeloma_all$outcome == "myeloma-finngen" ,]
mr_olink_protein_finngen_myeloma <-  mr_olink_protein_finngen_myeloma[grep("UKB-EU", mr_olink_protein_finngen_myeloma$exposure), ]
mr_olink_protein_finngen_myeloma$Gene <- sub("^(.*?)_.*$", "\\1", mr_olink_protein_finngen_myeloma$exposure)
mr_olink_protein_finngen_myeloma$UniProt <- sub("^[^_]*_(.*?)_.*$", "\\1", mr_olink_protein_finngen_myeloma$exposure)
mr_olink_protein_finngen_myeloma$b_ci <- 1.96 * mr_olink_protein_finngen_myeloma$se
mr_olink_protein_finngen_myeloma$b_ci_upper <- mr_olink_protein_finngen_myeloma$b + mr_olink_protein_finngen_myeloma$b_ci
mr_olink_protein_finngen_myeloma$b_ci_lower <- mr_olink_protein_finngen_myeloma$b - mr_olink_protein_finngen_myeloma$b_ci
mr_olink_protein_finngen_myeloma$or <- exp(mr_olink_protein_finngen_myeloma$b)
mr_olink_protein_finngen_myeloma$or_lower_ci <- exp(mr_olink_protein_finngen_myeloma$b - (1.96 * mr_olink_protein_finngen_myeloma$se)) 
mr_olink_protein_finngen_myeloma$or_upper_ci <- exp(mr_olink_protein_finngen_myeloma$b + (1.96 * mr_olink_protein_finngen_myeloma$se))
mr_olink_protein_finngen_myeloma <- merge(
  harmonised_data, 
  mr_olink_protein_finngen_myeloma, 
  by.x = c("exposure", "outcome"), 
  by.y = c("exposure", "outcome"), 
  all.y = TRUE
) %>%
 distinct()
#write.xlsx(mr_olink_protein_finngen_myeloma, file = "tables/mr_olink_protein_myeloma_finngen.xlsx", rowNames = F, colNames = T)

## Read in file to link in protein names for deCODE somalogic panel
prot_names <- readxl::read_xlsx(path = "Protein_myeloma_MR/ferkingstad_supplement.xlsx")
colnames(prot_names) <- gsub("\\s", "", colnames(prot_names))
prot_names <- prot_names[!duplicated(prot_names$Gene), ]
prot_names <- prot_names[,c("SeqId", "Protein(shortname)", "Protein(fullname)", "Gene", "UniProt")]
mr_soma_protein_ukb_myeloma$SeqId <- sub("^(.*?)_(.*?)_.*$", "\\1_\\2", mr_soma_protein_ukb_myeloma$exposure)
mr_soma_protein_ukb_myeloma_full <- merge(mr_soma_protein_ukb_myeloma, prot_names, by = "SeqId", all.x = TRUE, all.y = FALSE)
extract_string <- function(x) {
  parts <- unlist(strsplit(as.character(x), "_"))
  if(length(parts) >= 3) {
    return(parts[3])
  } else {
    return(NA)
  }
}
mr_soma_protein_ukb_myeloma_full$Gene <- sapply(mr_soma_protein_ukb_myeloma_full$exposure, extract_string)
mr_soma_protein_ukb_myeloma_full$b_ci <- 1.96 * mr_soma_protein_ukb_myeloma_full$se
mr_soma_protein_ukb_myeloma_full$b_ci_upper <- mr_soma_protein_ukb_myeloma_full$b + mr_soma_protein_ukb_myeloma_full$b_ci
mr_soma_protein_ukb_myeloma_full$b_ci_lower <- mr_soma_protein_ukb_myeloma_full$b - mr_soma_protein_ukb_myeloma_full$b_ci
mr_soma_protein_ukb_myeloma_full$or <- exp(mr_soma_protein_ukb_myeloma_full$b)
mr_soma_protein_ukb_myeloma_full$or_lower_ci <- exp(mr_soma_protein_ukb_myeloma_full$b - (1.96 * mr_soma_protein_ukb_myeloma_full$se)) 
mr_soma_protein_ukb_myeloma_full$or_upper_ci <- exp(mr_soma_protein_ukb_myeloma_full$b + (1.96 * mr_soma_protein_ukb_myeloma_full$se))
mr_soma_protein_ukb_myeloma_full <- merge(
  harmonised_data, 
  mr_soma_protein_ukb_myeloma_full, 
  by.x = c("exposure", "outcome"), 
  by.y = c("exposure", "outcome"), 
  all.y = TRUE
) %>%
  distinct()
write.xlsx(mr_soma_protein_ukb_myeloma_full, file = "tables/mr_soma_protein_myeloma_meta.xlsx", rowNames = F, colNames = T)

## Plot of deCODE and myeloma (UKB-finngen
mr_soma_protein_ukb_myeloma_subset <- subset(mr_soma_protein_ukb_myeloma_full, pval < 0.05)
mr_soma_protein_ukb_myeloma_subset <- mr_soma_protein_ukb_myeloma_subset[order(mr_soma_protein_ukb_myeloma_subset$pval), ] ## 53
length(which(mr_soma_protein_ukb_myeloma_subset$b > 0)) ## 20
length(which(mr_soma_protein_ukb_myeloma_subset$b < 0)) ## 33

mr_soma_protein_ukb_myeloma_subset$Gene <- with(mr_soma_protein_ukb_myeloma_subset,
                                                                  reorder(Gene, pval, FUN = function(x) -log10(x)))

#pdf(file = "figures/forest_mr_soma_ukbfinngen_myeloma.pdf", height = 8, width = 10)
ggplot(mr_soma_protein_ukb_myeloma_subset, aes(x = b, xmin = b_ci_lower, xmax = b_ci_upper, y = Gene)) +
  geom_point() +
  geom_errorbarh(aes(xmin = b_ci_lower, xmax = b_ci_upper), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Log odds of MM ± 95% CI per normalized SD unit higher protein", y = "Protein") +
  coord_cartesian(xlim = c(-11, NA)) +  # Set the lower limit of the x-axis to -11
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, angle = 0, hjust = 1),  
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1), 
    axis.title.x = element_text(size = 16),                       
    axis.title.y = element_text(size = 16)                        
  )

dev.off()

## Plot of Olink UKB and myeloma FinnGen
mr_olink_protein_finngen_myeloma_subset <- subset(mr_olink_protein_finngen_myeloma, pval < 0.05) ## 78
mr_olink_protein_finngen_myeloma_subset <- mr_olink_protein_finngen_myeloma_subset[order(mr_olink_protein_finngen_myeloma_subset$pval), ]
length(which(mr_olink_protein_finngen_myeloma_subset$b > 0)) ## 39
length(which(mr_olink_protein_finngen_myeloma_subset$b < 0)) ## 39

# Read in olink file to give full name
olink_names <- readxl::read_excel("003_protein_myeloma/olink-explore-3072-assay-list-2023-06-08.xlsx")
mr_olink_protein_finngen_myeloma_subset_full <- merge(mr_olink_protein_finngen_myeloma_subset, olink_names, by.x = "Gene", by.y = "Gene name" )
mr_olink_protein_finngen_myeloma_subset_full$Gene <- with(mr_olink_protein_finngen_myeloma_subset_full,
                                                                  reorder(Gene, pval, FUN = function(x) -log10(x)))

#pdf(file = "figures/forest_mr_olink_finngen_myeloma.pdf", height = 10, width = 10)
ggplot(mr_olink_protein_finngen_myeloma_subset_full, aes(x = b, xmin = b_ci_lower, xmax = b_ci_upper, y = Gene)) +
  geom_point() +
  geom_errorbarh(aes(xmin = b_ci_lower, xmax = b_ci_upper), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Log odds of MM risk ± 95% CI per normalised SD unit higher protein", y = "Protein") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, angle = 0, hjust = 1),  
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1), 
    axis.title.x = element_text(size = 16),                       
    axis.title.y = element_text(size = 16)                        
  )
dev.off()

## Olink UKB and myeloma UKB/FinnGen
mr_olink_protein_finngen_ukb_myeloma <-  mr_olink_protein_myeloma_all[mr_olink_protein_myeloma_all$outcome == "myeloma-UKB-finngen" ,]
mr_olink_protein_finngen_ukb_myeloma <-  mr_olink_protein_finngen_ukb_myeloma[grep("UKB-EU", mr_olink_protein_finngen_ukb_myeloma$exposure), ]
mr_olink_protein_finngen_ukb_myeloma$Gene <- sub("^(.*?)_.*$", "\\1", mr_olink_protein_finngen_ukb_myeloma$exposure)
mr_olink_protein_finngen_ukb_myeloma$UniProt <- sub("^[^_]*_(.*?)_.*$", "\\1", mr_olink_protein_finngen_ukb_myeloma$exposure)
mr_olink_protein_finngen_ukb_myeloma$b_ci <- 1.96 * mr_olink_protein_finngen_ukb_myeloma$se
mr_olink_protein_finngen_ukb_myeloma$b_ci_upper <- mr_olink_protein_finngen_ukb_myeloma$b + mr_olink_protein_finngen_ukb_myeloma$b_ci
mr_olink_protein_finngen_ukb_myeloma$b_ci_lower <- mr_olink_protein_finngen_ukb_myeloma$b - mr_olink_protein_finngen_ukb_myeloma$b_ci
mr_olink_protein_finngen_ukb_myeloma$or <- exp(mr_olink_protein_finngen_ukb_myeloma$b)
mr_olink_protein_finngen_ukb_myeloma$or_lower_ci <- exp(mr_olink_protein_finngen_ukb_myeloma$b - (1.96 * mr_olink_protein_finngen_ukb_myeloma$se)) 
mr_olink_protein_finngen_ukb_myeloma$or_upper_ci <- exp(mr_olink_protein_finngen_ukb_myeloma$b + (1.96 * mr_olink_protein_finngen_ukb_myeloma$se))
mr_olink_protein_finngen_ukb_myeloma <- merge(
  harmonised_data, 
  mr_olink_protein_finngen_ukb_myeloma, 
  by.x = c("exposure", "outcome"), 
  by.y = c("exposure", "outcome"), 
  all.y = TRUE
) %>%
  distinct()
#write.xlsx(mr_olink_protein_finngen_ukb_myeloma, file = "tables/mr_olink_protein_finngen_ukb_meta_myeloma.xlsx", rowNames = F, colNames = T)

## Scatter plot of the results using both myeloma GWAS with UKB olink protein exposure
sensitivity <- merge(mr_olink_protein_finngen_myeloma_subset, mr_olink_protein_finngen_ukb_myeloma, by = "Gene")

#pdf("figures/scatter_olink_myeloma_outcomes_comparison.pdf", height = 6, width = 6)
gg <- ggplot(sensitivity, aes(x = b.x, y = b.y)) +
  geom_point() +
  geom_errorbar(aes(x = b.x, ymin = b.y - se.y, ymax = b.y + se.y), width = 0.1, alpha = 0.3) +
  geom_errorbarh(aes(y = b.y, xmin = b.x - se.x, xmax = b.x + se.x), height = 0.1, alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(x = "UKB Olink protein and FinnGen myeloma GWAS",
       y = "UKB Olink protein and FinnGen/UKB myeloma GWAS") +
  theme_minimal()
dev.off() ## pearson correlation coefficient is 0.77

## Are all the proteins where p<0.05 using Olink and FinnGen also p<0.05 using the meta-analysis of MM GWAS?
mr_olink_protein_finngen_ukb_myeloma_subset <- subset(mr_olink_protein_finngen_ukb_myeloma, pval < 0.05)
w <- which(mr_olink_protein_finngen_myeloma_subset$Gene %in% mr_olink_protein_finngen_ukb_myeloma_subset$Gene)

## Are effect estimates in the same direction?
length(which(sensitivity$b.x > 0 & sensitivity$b.y > 0)) ## 653
length(which(sensitivity$b.x < 0 & sensitivity$b.y < 0)) ## 582

## What is the overlap in the proteins in both MR studies?
w <- which(mr_soma_protein_ukb_myeloma_full$Gene %in% mr_olink_protein_finngen_myeloma$Gene)
mr_olink_protein_finngen_myeloma <-  mr_olink_protein_myeloma_all[mr_olink_protein_myeloma_all$outcome == "myeloma-finngen" ,]
mr_olink_protein_finngen_myeloma <-  mr_olink_protein_finngen_myeloma[grep("UKB-EU", mr_olink_protein_finngen_myeloma$exposure), ]
mr_olink_protein_finngen_myeloma$Gene <- sub("^(.*?)_.*$", "\\1", mr_olink_protein_finngen_myeloma$exposure)
mr_olink_protein_finngen_myeloma$UniProt <- sub("^[^_]*_(.*?)_.*$", "\\1", mr_olink_protein_finngen_myeloma$exposure)
mr_olink_protein_finngen_myeloma$b_ci <- 1.96 * mr_olink_protein_finngen_myeloma$se
mr_olink_protein_finngen_myeloma$b_ci_upper <- mr_olink_protein_finngen_myeloma$b + mr_olink_protein_finngen_myeloma$b_ci
mr_olink_protein_finngen_myeloma$b_ci_lower <- mr_olink_protein_finngen_myeloma$b - mr_olink_protein_finngen_myeloma$b_ci
mr_olink_protein_finngen_myeloma$or <- exp(mr_olink_protein_finngen_myeloma$b)
mr_olink_protein_finngen_myeloma$or_lower_ci <- exp(mr_olink_protein_finngen_myeloma$b - (1.96 * mr_olink_protein_finngen_myeloma$se)) 
mr_olink_protein_finngen_myeloma$or_upper_ci <- exp(mr_olink_protein_finngen_myeloma$b + (1.96 * mr_olink_protein_finngen_myeloma$se))
w2 <- which(mr_olink_protein_finngen_myeloma$Gene %in% mr_soma_protein_ukb_myeloma_full$Gene)
length(w) ## 441
overlap_gene_names <- mr_soma_protein_ukb_myeloma_full$Gene[w]

## Do these have the same instruments?
exposure_dat <- read.table("003_protein_myeloma/data_exposure.txt", header = T)
## Exploring exposure data
exposure_dat_UKB_EU <-  exposure_dat[grep("UKB-EU", exposure_dat$exposure), ] ## 1860 proteins, median F stat 742
exposure_dat_UKB_EU <- exposure_dat_UKB_EU  %>%
  mutate(Gene = sub("^(.*?)_.*", "\\1", exposure))
#write.xlsx(exposure_dat_UKB_EU, file = "tables/exposure_data_ukb_olink.xlsx", colNames = T, rowNames = F)

exposure_dat_deCODE <-  exposure_dat[grep("DECODE", exposure_dat$exposure), ] ## 1192 proteins, median F stat 321
exposure_dat_deCODE <- exposure_dat_deCODE  %>%
  mutate(Gene = sub("^[^_]+_[^_]+_([^_]+).*", "\\1", exposure))
#write.xlsx(exposure_dat_deCODE, file = "tables/exposure_data_decode_somalogic.xlsx", colNames = T, rowNames = F)

rsid_dataframe_rows <- which(exposure_dat_UKB_EU$Gene %in% overlap_gene_names)
rsid_dataframe <- exposure_dat_UKB_EU[rsid_dataframe_rows,]

rsid_dataframe_rows_2 <- which(exposure_dat_deCODE$Gene %in% overlap_gene_names)
rsid_dataframe_2 <- exposure_dat_deCODE[rsid_dataframe_rows_2,]

## merge rsid dataframes
rsid <- merge(x = rsid_dataframe, y = rsid_dataframe_2, by = "Gene")
rsid$same_snp <- ifelse(rsid$SNP.x == rsid$SNP.y, "Yes", "No")
table(rsid$same_snp) # yes = 157, no = 285
rsid$same_chrom <- ifelse(rsid$chr.exposure.x == rsid$chr.exposure.y, "Yes", "No")
table(rsid$same_chrom) # yes = 442

prot_both_soma <- mr_soma_protein_ukb_myeloma_full[w,]
prot_both_olink <- mr_olink_protein_finngen_myeloma[w2,]

## Scatter plot of the comparison of these proteins
mr_comparison <- merge(prot_both_olink, prot_both_soma, by = "Gene", all.x = T, all.y = F) 
#write.xlsx(mr_comparison, file = 'tables/mr_comparison_discovery_replication.xlsx', rowNames = F, colNames = T)

duplicated_rows <- duplicated(mr_comparison$Gene)
mr_comparison <- mr_comparison[!duplicated_rows, ] ## 440

gg <- ggplot(mr_comparison, aes(x = b.x, y = b.y)) +
  geom_point() +
  geom_errorbar(aes(x = b.x, ymin = b.y - se.y, ymax = b.y + se.y), width = 0.1, alpha = 0.3) +
  geom_errorbarh(aes(y = b.y, xmin = b.x - se.x, xmax = b.x + se.x), height = 0.1, alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(x = "Discovery protein-multiple myeloma MR estimate",
       y = "Replication protein-multiple myeloma MR estimate") +
  theme_minimal()

# Calculate correlation coefficient and p-value
correlation <- cor.test(mr_comparison$b.x, mr_comparison$b.y)
r <- correlation$estimate
p_value <- correlation$p.value

# Format p-value
p_value_text <- ifelse(p_value < 0.001, "< 0.001", sprintf("%.3f", p_value))

# Add text annotation
#pdf(file = "figures/scatter_protein_MM_MR_estimates.pdf", height = 8, width = 8)
gg + annotate("text", x = max(mr_comparison$b.x), y = max(mr_comparison$b.y),
              label = paste("r =", round(r, 2), ", p =", p_value_text), 
              hjust = 0, vjust = 1)
dev.off()

positives_both <- which(mr_comparison$b.x > 0 & mr_comparison$b.y > 0) ## 167
negatives_both <- which(mr_comparison$b.x < 0 & mr_comparison$b.y < 0) ## 135

## Do any proteins have MR evidence in both MR studies
mr_soma_protein_ukb_myeloma_subset$Gene <- as.character(mr_soma_protein_ukb_myeloma_subset$Gene)
w <- which(mr_soma_protein_ukb_myeloma_subset$Gene %in% mr_olink_protein_finngen_myeloma_subset_full$Gene)
mr_soma_protein_ukb_myeloma_subset$Gene[w]
## [1] DPT    CRYBB1 OBP2B  GCLM   IL18BP KDR    CRYGD 

## Plot the comparison of estimates
common_prots <- mr_soma_protein_ukb_myeloma_subset$Gene[w]
common_prot_soma <- subset(mr_soma_protein_ukb_myeloma_subset, Gene %in% common_prots)
common_prot_soma$Study <- "SomaLogic protein GWAS, FinnGen + UKB myeloma meta-analysis of GWAS"
common_prot_olink <- subset(mr_olink_protein_finngen_myeloma_subset, Gene %in% common_prots)
common_prot_olink$Study <- "Olink protein GWAS, FinnGen myeloma GWAS"

## Restrict to column columns and rbind
common_columns <- intersect(names(common_prot_soma), names(common_prot_olink))
restricted_soma <- common_prot_soma[, common_columns]
restricted_olink <- common_prot_olink[, common_columns]

# Combine the two dataframes vertically (rbind)
combined_dataframe_prot <- rbind(restricted_soma, restricted_olink)

## Plot both MR estimates for proteins with MR evidence for an effect
#pdf(file = "figures/forest_mr_consistent_prots_meta.pdf", height = 6, width = 8)
ggplot(combined_dataframe_prot, aes(x = b, xmin = b_ci_lower, xmax = b_ci_upper, y = Gene, color = Study)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = b_ci_lower, xmax = b_ci_upper), height = 0.1, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Log odds of MM risk ± 95% CI per normalized SD unit higher protein", y = "Protein") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),  # Rotate the y-axis labels
    axis.text.x = element_text(angle = 45, hjust = 1),   # Rotate the x-axis labels
    legend.position = "top"  # Move the legend to the top
  ) +
  guides(color = guide_legend(nrow = 2, title = "Study", label.theme = element_text(lineheight = 0.9))) +
  scale_color_discrete(labels = function(labels) sapply(labels, function(label) paste(strwrap(label, width = 20), collapse = "\n")))
dev.off()

## Myeloma -> protein results
mr_myeloma_protein <- read.delim("004_myeloma_protein/data_mr.txt")
mr_myeloma_protein$or <- exp(mr_myeloma_protein$b)
mr_myeloma_protein$or_lower_ci <- exp(mr_myeloma_protein$b - (1.96 * mr_myeloma_protein$se)) 
mr_myeloma_protein$or_upper_ci <- exp(mr_myeloma_protein$b + (1.96 * mr_myeloma_protein$se))

## Restrict to ukb/finngen myeloma deCODE proteins
mr_ukb_myeloma <- mr_myeloma_protein[mr_myeloma_protein$exposure == "myeloma-UKB-finngen" ,]
mr_ukb_myeloma_soma_protein <- mr_ukb_myeloma[grepl("^\\d", mr_ukb_myeloma[["id.outcome"]]), ]

## Merge in protein names
mr_ukb_myeloma_soma_protein$SeqId <- sub("^(.*?)_(.*?)_.*$", "\\1_\\2", mr_ukb_myeloma_soma_protein$outcome)
mr_ukb_myeloma_soma_protein_full <- merge(mr_ukb_myeloma_soma_protein, prot_names, by = "SeqId", all.x = TRUE, all.y = FALSE)
mr_ukb_myeloma_soma_protein_full$b_ci <- 1.96 * mr_ukb_myeloma_soma_protein_full$se
mr_ukb_myeloma_soma_protein_full$b_ci_upper <- mr_ukb_myeloma_soma_protein_full$b + mr_ukb_myeloma_soma_protein_full$b_ci
mr_ukb_myeloma_soma_protein_full$b_ci_lower <- mr_ukb_myeloma_soma_protein_full$b - mr_ukb_myeloma_soma_protein_full$b_ci
mr_ukb_myeloma_soma_protein_full$Gene <- with(mr_ukb_myeloma_soma_protein_full,  reorder(Gene, pval, FUN = function(x) -log10(x)))
#write.xlsx(mr_ukb_myeloma_soma_protein_full, file = "tables/mr_myeloma_meta_soma.xlsx", rowNames = F, colNames = T)

## Restrict to those where p<0.05
mr_ukb_myeloma_soma_protein_full_subset <- subset(mr_ukb_myeloma_soma_protein_full, pval < 0.05)
mr_ukb_myeloma_soma_protein_full_subset <- mr_ukb_myeloma_soma_protein_full_subset[order(mr_ukb_myeloma_soma_protein_full_subset$pval), ] ## 44
length(which(mr_ukb_myeloma_soma_protein_full_subset$b > 0)) ## 35
length(which(mr_ukb_myeloma_soma_protein_full_subset$b < 0)) ## 9

## FinnGen myeloma and UKB olink results
mr_myeloma_all_olink_protein <- mr_myeloma_protein[str_detect(mr_myeloma_protein$id.outcome, "^[A-Z]" ), ]
mr_finngen_myeloma_olink_protein <-  mr_myeloma_all_olink_protein[mr_myeloma_all_olink_protein$exposure == "myeloma-finngen" ,]
mr_finngen_myeloma_olink_protein <-  mr_finngen_myeloma_olink_protein[grep("UKB-EU", mr_finngen_myeloma_olink_protein$outcome), ]
mr_finngen_myeloma_olink_protein$Gene <- sub("^(.*?)_.*$", "\\1", mr_finngen_myeloma_olink_protein$outcome)
mr_finngen_myeloma_olink_protein$UniProt <- sub("^[^_]*_(.*?)_.*$", "\\1", mr_finngen_myeloma_olink_protein$outcome)
mr_finngen_myeloma_olink_protein$b_ci <- 1.96 * mr_finngen_myeloma_olink_protein$se
mr_finngen_myeloma_olink_protein$b_ci_upper <- mr_finngen_myeloma_olink_protein$b + mr_finngen_myeloma_olink_protein$b_ci
mr_finngen_myeloma_olink_protein$b_ci_lower <- mr_finngen_myeloma_olink_protein$b - mr_finngen_myeloma_olink_protein$b_ci
#write.xlsx(mr_finngen_myeloma_olink_protein, file = "tables/mr_myeloma_finngen_olink.xlsx", rowNames = F, colNames = T)

mr_finngen_myeloma_olink_protein_subset <- subset(mr_finngen_myeloma_olink_protein, pval < 0.05)
mr_finngen_myeloma_olink_protein_subset <- mr_finngen_myeloma_olink_protein_subset[order(mr_finngen_myeloma_olink_protein_subset$pval), ] ## 136
length(unique(mr_finngen_myeloma_olink_protein_subset$outcome)) ## 136
length(which(mr_finngen_myeloma_olink_protein_subset$b > 0)) ## 67
length(which(mr_finngen_myeloma_olink_protein_subset$b < 0)) ## 69

## Which proteins are associated in both directions in each MR
## UKB/Finngen myeloma and deCODE somalogic proteins
both_directions_1 <- which(mr_soma_protein_ukb_myeloma_subset$Gene %in% mr_ukb_myeloma_soma_protein_full_subset$Gene)
print(mr_soma_protein_ukb_myeloma_subset$Gene[both_directions_1]) ## MMP9
r <- mr_soma_protein_ukb_myeloma_subset$Gene[both_directions_1]
w <- which(mr_soma_protein_ukb_myeloma_subset$Gene %in% r)

## Remove this protein from plot
mr_soma_protein_ukb_myeloma_subset_rr <- mr_soma_protein_ukb_myeloma_subset[-w,]

#pdf(file = "figures/forest_mr_soma_ukbfinngen_myeloma_reverse_removed.pdf", height = 8, width = 10)
ggplot(mr_soma_protein_ukb_myeloma_subset, aes(x = b, xmin = b_ci_lower, xmax = b_ci_upper, y = Gene)) +
  geom_point() +
  geom_errorbarh(aes(xmin = b_ci_lower, xmax = b_ci_upper), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Log odds of MM ± 95% CI per normalized SD unit higher protein", y = "Protein") +
  coord_cartesian(xlim = c(-11, NA)) +  # Set the lower limit of the x-axis to -11
  theme_minimal() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),  # Rotate the y-axis labels
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate the x-axis labels
  )
dev.off()
  
both_directions_2 <- which(mr_olink_protein_finngen_myeloma_subset$Gene %in% mr_finngen_myeloma_olink_protein_subset$Gene)
print(mr_olink_protein_finngen_myeloma_subset$Gene[both_directions_2]) ##  [1] "KLK15" 

## Repeat for the replication MR - remove proteins from plot that have evidence for an effect in the reverse MR
r2 <- as.character(mr_olink_protein_finngen_myeloma_subset_full$Gene[both_directions_2])
w2 <- which(mr_olink_protein_finngen_myeloma_subset_full$Gene %in% r2)

## Remove this protein from plot
mr_olink_protein_finngen_myeloma_subset_full_rr <- mr_olink_protein_finngen_myeloma_subset_full[-w2,] ## This removes one protein- OBP2B

#pdf(file = "figures/forest_mr_olink_finngen_myeloma_rr.pdf", height = 8, width = 10)
ggplot(mr_olink_protein_finngen_myeloma_subset_full_rr, aes(x = b, xmin = b_ci_lower, xmax = b_ci_upper, y = Gene)) +
  geom_point() +
  geom_errorbarh(aes(xmin = b_ci_lower, xmax = b_ci_upper), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Log odds of MM risk ± 95% CI per normalised SD unit higher protein", y = "Protein") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),  # Rotate the y-axis labels
    axis.text.x = element_text(angle = 45, hjust = 1)   # Rotate the x-axis labels
  )
dev.off()

## Combine the forward and reverse MR estimates for the 7 proteins with consistent evidence across both MR analyses
w <- which(mr_finngen_myeloma_olink_protein$Gene %in% common_prots)
mr_finngen_myeloma_olink_protein_common <- mr_finngen_myeloma_olink_protein[w,]
w <- which(mr_finngen_myeloma_olink_protein_common$method == "Inverse variance weighted")
mr_finngen_myeloma_olink_protein_common_ivw <- mr_finngen_myeloma_olink_protein_common[w,]

w2 <- which(mr_ukb_myeloma_soma_protein_full$Gene %in% common_prots)
mr_ukb_myeloma_soma_protein_full_common <- mr_ukb_myeloma_soma_protein_full[w2,]

## Combine these MR estimates with the forward MR estimates
## Restrict to column columns and rbind
common_columns <- intersect(names(combined_dataframe_prot), names(mr_finngen_myeloma_olink_protein_common_ivw))
restricted_olink_reverse <- mr_finngen_myeloma_olink_protein_common_ivw[, common_columns]
restricted_olink_reverse$Study <- "Olink protein GWAS, FinnGen myeloma GWAS"
restricted_olink_reverse$Direction <- "Reverse"
restricted_somalogic_reverse <- mr_ukb_myeloma_soma_protein_full_common[, common_columns]
restricted_somalogic_reverse$Study <- "SomaLogic protein GWAS, FinnGen + UKB myeloma meta-analysis of GWAS"
restricted_somalogic_reverse$Direction <- "Reverse"
combined_dataframe_prot2 <- combined_dataframe_prot[, common_columns]
combined_dataframe_prot2$Study[1:length(common_prots)] <- "SomaLogic protein GWAS, FinnGen + UKB myeloma meta-analysis of GWAS"
combined_dataframe_prot2$Study[(length(common_prots)+1):(length(common_prots)+length(common_prots))] <- "Olink protein GWAS, FinnGen myeloma GWAS"
combined_dataframe_prot2$Direction <- "Forward"

# Combine the two dataframes vertically (rbind)
combined_dataframe_prot2 <- rbind(combined_dataframe_prot2, restricted_olink_reverse, restricted_somalogic_reverse)
combined_dataframe_prot2$b_ci_upper <- combined_dataframe_prot2$b + 1.96 * combined_dataframe_prot2$se
combined_dataframe_prot2$b_ci_lower <- combined_dataframe_prot2$b - 1.96 * combined_dataframe_prot2$se

# Plot these results - those that are consistent, both forward and reverse MR
#pdf(file = "figures/forest_mr_consistent_prots_meta_forward_reverse.pdf", height = 6, width = 8)
ggplot(combined_dataframe_prot2, aes(x = b, xmin = b_ci_lower, xmax = b_ci_upper, y = Gene, color = Direction, shape = Study)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = b_ci_lower, xmax = b_ci_upper), height = 0.1, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Log odds of MM risk ± 95% CI per normalized SD unit higher protein(forward)\n or difference in protein levels (SDs + 95% CI) per effect allele increase in MM risk (reverse)", y = "Protein") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),  # Rotate the y-axis labels
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate the x-axis labels
    legend.position = "top"  # Move the legend to the top
  ) +
  guides(
    color = guide_legend(nrow = 2, title = "Direction", label.theme = element_text(lineheight = 0.9)),
    shape = guide_legend(nrow = 2, title = "Study", label.theme = element_text(lineheight = 0.9))
  ) +
  scale_color_discrete(labels = function(labels) sapply(labels, function(label) paste(strwrap(label, width = 20), collapse = "\n")))

dev.off()

## Steiger filtering - forward MR
steiger_prot_myeloma <- read.delim("003_protein_myeloma/data_steiger.txt")

## Steiger protein and myeloma direction, deCODE somalogic protien, ukb+finngen myeloma
steiger_soma_protein_ukb_myeloma <- steiger_prot_myeloma[steiger_prot_myeloma$outcome == "myeloma-UKB-finngen" ,]
steiger_soma_protein_ukb_myeloma <- subset(steiger_soma_protein_ukb_myeloma, grepl("DECODE", exposure)) 
#write.xlsx(steiger_soma_protein_ukb_myeloma, file = "tables/steiger_soma_protein_ukbfinngen_myeloma.xlsx", colNames = T, rowNames = F)

## How many have a FALSE for the forward direction
false_steiger_1 <- which(steiger_soma_protein_ukb_myeloma$correct_causal_direction == "FALSE")
steiger_soma_protein_ukb_myeloma$exposure[false_steiger_1] # none - all true

## Steiger protein and myeloma direction, ukb protein, finngen myeloma
steiger_olink_protein_myeloma_all <- steiger_prot_myeloma[str_detect(steiger_prot_myeloma$exposure, "^[A-Z]" ), ]
steiger_olink_protein_finngen_myeloma <-  steiger_olink_protein_myeloma_all[steiger_olink_protein_myeloma_all$outcome == "myeloma-finngen" ,]
steiger_olink_protein_finngen_myeloma <-  steiger_olink_protein_finngen_myeloma[grep("UKB-EU", steiger_olink_protein_finngen_myeloma$exposure), ]
steiger_olink_protein_finngen_myeloma$Gene <- sub("^(.*?)_.*$", "\\1", steiger_olink_protein_finngen_myeloma$exposure)
steiger_olink_protein_finngen_myeloma$UniProt <- sub("^[^_]*_(.*?)_.*$", "\\1", steiger_olink_protein_finngen_myeloma$exposure)
#write.xlsx(steiger_olink_protein_finngen_myeloma, file = "tables/steiger_olink_protein_finngen_myeloma.xlsx", colNames = T, rowNames = F)

## How many have a FALSE for the forward direction
false_steiger_2 <- which(steiger_olink_protein_finngen_myeloma$correct_causal_direction == "FALSE")
steiger_olink_protein_finngen_myeloma$exposure[false_steiger_2] # None
table(steiger_olink_protein_finngen_myeloma$correct_causal_direction) ## All 1622 are true.

## Steiger filtering - reverse MR 
steiger_myeloma_prot <- read.delim("004_myeloma_protein/data_steiger.txt")

## uk myeloma and deCODE somalogic
steiger_ukb_myeloma <- steiger_myeloma_prot[steiger_myeloma_prot$exposure == "myeloma-UKB-finngen" ,]
steiger_ukb_myeloma_soma_protein <- steiger_ukb_myeloma[grepl("^\\d", steiger_ukb_myeloma[["id.outcome"]]), ]
#write.xlsx(steiger_ukb_myeloma_soma_protein, file = "tables/steiger_ukb_myeloma_soma_protein.xlsx", colNames = T, rowNames = F)

## How many have a false for the reverse direction
false_steiger_3 <- which(steiger_ukb_myeloma_soma_protein$correct_causal_direction == "FALSE")
steiger_ukb_myeloma_soma_protein$exposure[false_steiger_3] # None
table(steiger_ukb_myeloma_soma_protein$correct_causal_direction) ## 4907 are TRUE

## finngen myeloma and ukbiobank protein
steiger_myeloma_all_olink_protein <- steiger_myeloma_prot[str_detect(steiger_myeloma_prot$outcome, "^[A-Z]" ), ]
steiger_finngen_myeloma_olink_protein <-  steiger_myeloma_all_olink_protein[steiger_myeloma_all_olink_protein$exposure == "myeloma-finngen" ,]
steiger_finngen_myeloma_olink_protein <-  steiger_finngen_myeloma_olink_protein[grep("UKB-EU", steiger_finngen_myeloma_olink_protein$id.outcome), ]
steiger_finngen_myeloma_olink_protein$Gene <- sub("^(.*?)_.*$", "\\1", steiger_finngen_myeloma_olink_protein$outcome)
steiger_finngen_myeloma_olink_protein$UniProt <- sub("^[^_]*_(.*?)_.*$", "\\1", steiger_finngen_myeloma_olink_protein$outcome)

## How many have a FALSE for the forward direction
false_steiger_4 <- which(steiger_finngen_myeloma_olink_protein$correct_causal_direction == "FALSE")
steiger_finngen_myeloma_olink_protein$exposure[false_steiger_4] # None
table(steiger_finngen_myeloma_olink_protein$correct_causal_direction) ## All 2930 are true.

## How many are true in both directions??
steiger_both_ways <- which(steiger_olink_protein_finngen_myeloma$Gene %in% steiger_finngen_myeloma_olink_protein$Gene)
list <- steiger_olink_protein_finngen_myeloma$Gene[steiger_both_ways] 
length(list) ## 1616

## Of these how many were associated in forward MR? 
steiger_mr_consistent_forward <- which(mr_olink_protein_finngen_myeloma_subset_full$Gene %in% list)
length(steiger_mr_consistent_forward)  ## 78
tab <- mr_olink_protein_finngen_myeloma_subset_full[steiger_mr_consistent_forward,]

## Of these, how many were associated in reverse MR?
steiger_mr_consistent_reverse <- which(mr_finngen_myeloma_olink_protein_subset$Gene %in% list) 
length(steiger_mr_consistent_reverse) ## 82
tab2 <- mr_finngen_myeloma_olink_protein_subset[steiger_mr_consistent_reverse,]

## what is the overlap in these proteins
steiger_MR_reverse <- which(tab$Gene %in% tab2$Gene) 
tab$Gene[steiger_MR_reverse] ## 78

## Coloc results
coloc <- read.delim("coloc/001_results-coloc.txt", header = T)

## Split coloc results to each study and save
w <- grepl("DECODE;myeloma-UKB-finngen", coloc$id)
coloc_decode_ukb_myeloma_all <- coloc[w,]
#write.xlsx(coloc_decode_ukb_myeloma_all, file = "tables/coloc_meta_myeloma_soma_protein.xlsx", colNames = T, rowNames = F)

w2 <- grepl('european;myeloma-finngen', coloc$id)
coloc_olink_finngen_myeloma_all <- coloc[w2,]
#write.xlsx(coloc_olink_finngen_myeloma_all, file = "tables/coloc_olink_finngen_myeloma_all.xlsx", colNames = T, rowNames = F)

## Those with evidence of colocalisation
coloc_h4 <- subset(coloc, h4 > 0.8) 

## Split out into the separate MR results of interest
## decode and ieu ukb myeloma MR
w <- grepl("DECODE;myeloma-UKB-finngen", coloc_h4$id)
coloc_decode_ukb_myeloma <- coloc_h4[w,] ## KIAA1161

## Olink and finngen myeloma MR
w2 <- grepl('european;myeloma-finngen', coloc_h4$id)
coloc_olink_finngen_myeloma <- coloc_h4[w2,] ## none

## Restrict coloc results to the 7 proteins of interest
filtered_coloc <- coloc[grepl(paste(common_prots, collapse="|"), coloc$exposure), ]
filtered_coloc <- filtered_coloc[grepl("DECODE;myeloma-UKB-finngen", filtered_coloc$id), ]

## Restrict coloc results to the 7 proteins of interest - replication
filtered_coloc2 <- coloc[grepl(paste(common_prots, collapse="|"), coloc$exposure), ]
filtered_coloc2 <- filtered_coloc2[grepl('european;myeloma-finngen', filtered_coloc2$id), ]

## Compare results to published study by Wang et al.
wang <- read.xlsx("Protein_myeloma_MR/001_protein_myeloma/Wang_MR_results.xlsx")
decode_interval <- merge(wang, mr_soma_protein_ukb_myeloma_full, by.x = "Exposures", by.y = "Gene", all.x = T)
plot(decode_interval$OR_IVW, decode_interval$or)
