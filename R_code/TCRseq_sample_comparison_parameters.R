library(readr)
library(readxl)
library(dplyr)

# TCRseq-sample comparison matrix
## Run whole repertoire comparisons for INT TCR data

clone_data <- readRDS("~/Documents/BFX_proj/TCRseq_pipeline/_output/clone_data.rds")

comparison_matrix <- data.frame(sample_a = "P10-PB1.tsv",
                                sample_b = unique(clone_data$file),
                                id = "nucleic_acid",
                                d_a_method = "binomial")

differential_clone_abundance_results_folder <- "~/Documents/BFX_proj/TCRseq_pipeline/_output/differential_clone_abundance/"
 

