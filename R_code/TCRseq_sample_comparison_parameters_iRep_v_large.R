library(readr)
library(readxl)
library(dplyr)

# v_large
## individual file data exists in a folder that must be referenced

# output folder
differential_clone_abundance_results_folder <- "~/Documents/BFX_proj/TCRseq_pipeline/_output/differential_clone_abundance/"

# TCRseq-sample comparison matrix

### TCRseq-sample comparison matrix --------------------------------------------
manifest <- read_csv("~/Documents/BFX_proj/_data_private/INT_PCV/TCR/P101_iRep_summary/Moderna_mRNA4157P101_Q170423_Rao_sample_metadata.csv") # load sample manifest
colnames(manifest) <- gsub(" ", "_", colnames(manifest)) # clean up column names
patients_w_multiple_samples <- manifest$Screen_ID[duplicated(manifest$Screen_ID)] # identify patients with paired samples
manifest <- manifest %>% filter(Screen_ID %in% patients_w_multiple_samples) # 35 samples

comparison_matrix <- data.frame() # initialize comparison matrix
for(i in unique(manifest$Screen_ID)){ # for each subject:
  comp_mat_ <- data.frame() # initialize comparison matrix for a given subject
  
  man_ <- manifest[manifest$Screen_ID == i, ] # subset manifest to just that subject
  earliest_ <- man_[man_$LBDAT == min(man_$LBDAT), ] # identify earliest sample to use as reference for comparisons
  samp_a_ <- earliest_$Data_Name_for_TCRseq
  samp_a_ <- gsub(" ", "_", do.call(paste, expand.grid(samp_a_, c("hTRA", "hTRB", "hTRD", "hTRG"))))
  samp_a_ <- rep(samp_a_, each = 2)
  samp_a_ <- paste0(samp_a_, "_csv.rds")
  bc_a_ <- earliest_$Sample_Barcode
  bc_a_ <- gsub(" ", "_", do.call(paste, expand.grid(bc_a_, c("hTRA", "hTRB", "hTRD", "hTRG"))))
  bc_a_ <- rep(bc_a_, each = 2)
  
  for(j in 1:nrow(man_)){
    samp_b_ <- man_[j, ]$Data_Name_for_TCRseq
    samp_b_ <- gsub(" ", "_", do.call(paste, expand.grid(samp_b_, c("hTRA", "r_hTRA", "hTRB", "r_hTRB", "hTRD", "r_hTRD", "hTRG", "r_hTRG"))))
    samp_b_ <- paste0(samp_b_, "_csv.rds")
    bc_b_ <- man_[j, ]$Sample_Barcode
    bc_b_ <- gsub(" ", "_", do.call(paste, expand.grid(bc_b_, c("hTRA", "r_hTRA", "hTRB", "r_hTRB", "hTRD", "r_hTRD", "hTRG", "r_hTRG"))))
    
    comp_mat_ <- rbind(comp_mat_, data.frame(sample_a = samp_a_, sample_b = samp_b_, barcode_a = bc_a_, barcode_b = bc_b_))
  }
  comparison_matrix <- rbind(comparison_matrix, comp_mat_)
  rm(i, j, comp_mat_, man_, earliest_, samp_a_, samp_b_)
}

comparison_matrix$id <- "cdr3_v_j"
comparison_matrix$d_a_method <- "binomial"
comparison_matrix$clone_data_folder <- "temp"

rm(manifest, patients_w_multiple_samples, sample_barcode, sample_name)

saveRDS(comparison_matrix, paste0("~/Documents/BFX_proj/INT/processed_data/TCRseq/iRep_2023_12/comparison_matrix_", gsub("-", "_", Sys.Date()), ".rds"))





