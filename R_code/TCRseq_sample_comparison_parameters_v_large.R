library(readr)
library(readxl)
library(dplyr)

# v_large
## individual file data exists in a folder that must be referenced

# output folder
differential_clone_abundance_results_folder <- "~/Documents/BFX_proj/TCRseq_pipeline/_output/differential_clone_abundance/"

# TCRseq-sample comparison matrix
## Run whole repertoire comparisons for INT TCR data
clone_data <- readRDS("~/Documents/BFX_proj/TCRseq_pipeline/_output/clone_data.rds")

comparison_matrix <- data.frame(clone_data_folder = "~/Documents/BFX_proj/TCRseq_pipeline/_output/temp/condensed/",
                                sample_a = "P10-PB1_tsv",
                                sample_b = gsub("\\.", "_", unique(clone_data$file)),
                                id = "nucleic_acid",
                                d_a_method = "binomial")


 

