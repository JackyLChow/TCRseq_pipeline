library(readr)
library(readxl)
library(dplyr)

# comparison matrix

comparison_matrix <- rbind(data.frame(name_a = "P10-PB1.rds",
                                      name_b = paste0("P10-PB", 1:3, ".rds"),
                                      id = "nucleic_acid",
                                      d_a_method = "binomial"),
                           data.frame(name_a = "P9-PB1.rds",
                                      name_b = paste0("P9-PB", 1:3, ".rds"),
                                      id = "nucleic_acid",
                                      d_a_method = "binomial"))
                           
comparison_matrix$sample_a <- paste0("~/Documents/BFX_proj/TCRseq_pipeline/_input/counts/", comparison_matrix$name_a)
comparison_matrix$sample_b <- paste0("~/Documents/BFX_proj/TCRseq_pipeline/_input/counts/", comparison_matrix$name_b)

# for(i in gsub(".rds", ".tsv", comparison_matrix$sample_b)){
#   j <- gsub("~/Documents/BFX_proj/TCRseq_pipeline/_input/counts/", "", i)
#   j <- gsub(".tsv", ".rds", j)
#   f_ <- read.table(i, sep = "\t", header = T)
#   f_$count <- f_$count..templates.reads.
#   f_$nucleic_acid <- f_$nucleotide
#   saveRDS(f_, paste0("~/Documents/BFX_proj/TCRseq_pipeline/_input/counts/", j))
# }
 

