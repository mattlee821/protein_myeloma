# Exploring the role of circulating proteins in multiple myeloma risk: a Mendelian randomization study 

In this work we investigated the assocation between circulating proteins and risk of multiple myeloma (MM) using Mendelian randomisation and genetic colocalisation. The manuscript is available: [10.21203/rs.3.rs-4800219](https://www.researchsquare.com/article/rs-4800219/).

## [`scripts/`](https://github.com/mattlee821/protein_myeloma/tree/main/scripts)
Code is separated into directories and files and is labelled for the order in which it was run. For example, `001_instruments/001_instruments-proteins.R` was run before `002_outcomes/001_instruments-cancer_outcome-proteins.R`. Where a label number is the same, e.g., `002_outcomes/001_instruments-cancer_outcome-proteins.R` and `002_outcomes/001_instruments-proteins_outcome-cancer.R`, there is no order preference. Each directory and script is named for its specific purpose and scripts have extensive commenting to facilitate re-use.

* `001_instruments/`: provides the code for identifying and extracting genetic variants associated with proteins (`001_instruments-proteins.R`) and MM (`002_instruments-cancer.R`)
* `002_outcomes/`: provides the code for extracting the instruments obtained from `001_instruments/` for either MM `001_instruments-proteins_outcome-cancer.R` or proteins `001_instruments-cancer_outcome-proteins/`
* `003_protein-myeloma/`: provides the code `001_protein-myeloma.R` and the associated slurm submission script `001_protein-myeloma.sh` to perform the MR analyses (including Steiger test) of protein-MM risk
* `004_myeloma-protein/`: provides the code `001_MR.R` and the associated slurm submission script `001_MR.sh` to perform the MR analyses (including Steiger test) of MM risk-protein
* `005_colocalization/`: provides the code to extract windows for proteins (`001_exposure-window/`), to extract those windows from MM (`002_outcome-window/`), to ensure the same SNPs are present in both windows (`003_harmonise-window/`), and to perform colocalisation (`004_coloc/`) for different window sizes (e.g., `004_coloc/1mb/`).
