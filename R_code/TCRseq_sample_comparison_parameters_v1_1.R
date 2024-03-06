library(readr)
library(readxl)
library(dplyr)

# comparison matrix

comparison_matrix <- rbind(data.frame(name_a = "P10-PB1.tsv",
                                      name_b = paste0("P10-PB", 1:3, ".tsv"),
                                      id = "nucleic_acid",
                                      d_a_method = "binomial"),
                           data.frame(name_a = "P9-PB1.tsv",
                                      name_b = paste0("P9-PB", 1:3, ".tsv"),
                                      id = "nucleic_acid",
                                      d_a_method = "binomial"))
                           
comparison_matrix$sample_a <- paste0("~/Documents/BFX_proj/TCRseq_pipeline/_input/counts/", comparison_matrix$name_a)
comparison_matrix$sample_b <- paste0("~/Documents/BFX_proj/TCRseq_pipeline/_input/counts/", comparison_matrix$name_b)

 

