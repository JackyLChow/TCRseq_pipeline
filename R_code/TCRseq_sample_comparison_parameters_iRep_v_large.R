library(readr)
library(readxl)
library(dplyr)

# v_large
## individual file data exists in a folder that must be referenced

# output folder
differential_clone_abundance_results_folder <- "/"

# TCRseq-sample comparison matrix

### TCRseq-sample comparison matrix --------------------------------------------
manifest <- read_csv("") # load sample manifest
colnames(manifest) <- gsub(" ", "_", colnames(manifest)) # clean up column names
patients_w_multiple_samples <- manifest$Screen_ID[duplicated(manifest$Screen_ID)] # identify patients with paired samples
manifest <- manifest %>% filter(Screen_ID %in% patients_w_multiple_samples) # 35 samples

comparison_matrix <- data.frame() # initialize comparison matrix
for(i in unique(manifest$Screen_ID)){ # for each subject:
  comp_mat_ <- data.frame() # initialize comparison matrix for a given subject
  
  man_ <- manifest[manifest$Screen_ID == i, ] # subset manifest to just that subject
  earliest_ <- man_[grepl("Screen", man_$Visit), ] # identify earliest sample to use as reference for comparisons
  samp_a_ <- earliest_$Data_Name_for_TCRseq
  samp_a_ <- gsub(" ", "_", do.call(paste, expand.grid(samp_a_, paste0("PER_UMI_", c("hTRA", "hTRB", "hTRD", "hTRG")))))
  samp_a_ <- rep(samp_a_, each = 2)
  samp_a_ <- paste0(samp_a_, "_airr_csv.rds")
  bc_a_ <- earliest_$Sample_Barcode
  name_a_ <- gsub(" ", "_", do.call(paste, expand.grid(bc_a_, paste0("PER_UMI_", c("hTRA", "hTRB", "hTRD", "hTRG")))))
  name_a_ <- rep(name_a_, each = 2)
  
  for(j in 1:nrow(man_)){
    samp_b_ <- man_[j, ]$Data_Name_for_TCRseq
    samp_b_ <- gsub(" ", "_", do.call(paste, expand.grid(samp_b_, c("PER_UMI_hTRA", "r_PER_UMI_hTRA", "PER_UMI_hTRB", "r_PER_UMI_hTRB", "PER_UMI_hTRD", "r_PER_UMI_hTRD", "PER_UMI_hTRG", "r_PER_UMI_hTRG"))))
    samp_b_ <- paste0(samp_b_, "_airr_csv.rds")
    bc_b_ <- man_[j, ]$Sample_Barcode
    name_b_ <- gsub(" ", "_", do.call(paste, expand.grid(bc_b_, c("PER_UMI_hTRA", "r_PER_UMI_hTRA", "PER_UMI_hTRB", "r_PER_UMI_hTRB", "PER_UMI_hTRD", "r_PER_UMI_hTRD", "PER_UMI_hTRG", "r_PER_UMI_hTRG"))))
    
    comp_mat_ <- rbind(comp_mat_, data.frame(sample_a = samp_a_, sample_b = samp_b_, barcode_a = bc_a_, barcode_b = bc_b_, name_a = name_a_, name_b = name_b_))
  }
  comparison_matrix <- rbind(comparison_matrix, comp_mat_)
  rm(i, j, comp_mat_, man_, earliest_, samp_a_, samp_b_, bc_a_, bc_b_, name_a_, name_b_)
}

comparison_matrix$id <- "cdr3_v_j"
comparison_matrix$d_a_method <- "binomial"

head(comparison_matrix, 12)

### convert sample names to file paths -----------------------------------------
# get files paths for all condensed and downsampled data
files <- list.files("", full.names = T, recursive = T, pattern = "airr_csv.rds$")
files <- files[grepl("/downsampled/", files)]

# check sample_a 
## if sample_a non-replicate is missing, replace with replicate
## only one original sample missing: barcode 12N09536274
for(i in unique(comparison_matrix$sample_a)){
  if(sum(grepl(i, files)) != 1){
    cat(paste(i, "missing\n"))
    rep_ <- gsub("_TCR_", "_TCR_r_", i)
    cat(paste("replacing with", rep_, "\n"))
    comparison_matrix$sample_a[comparison_matrix$sample_a == i] <- rep_
    comparison_matrix$name_a[comparison_matrix$sample_a == i] <- gsub("PER_UMI", "r_PER_UMI", i)
    rm(rep_)
  }
  rm(i)
}

# check sample_b
## missing:
comparison_matrix$sample_b <- gsub("Not_Provided", "NotProvided", comparison_matrix$sample_b)
oof_sample_b <- c()
for(i in unique(comparison_matrix$sample_b)){
  if(sum(grepl(i, files)) != 1){
    cat(paste(i, "missing\n"))
    oof_sample_b  <- c(oof_sample_b , i)
  }
}
## filter out missing sample_b
comparison_matrix <- comparison_matrix[!comparison_matrix$sample_b %in% oof_sample_b, ]

# overwrite sample with file path
for(i in unique(comparison_matrix$sample_a)){
  path_ <- files[grepl(i, files)]
  comparison_matrix$sample_a[comparison_matrix$sample_a == i] <- path_
  rm(i, path_)
}
for(i in 1:nrow(comparison_matrix)){
  path_ <- files[grepl(comparison_matrix[i, "sample_b"], files)]
  comparison_matrix[i, "sample_b"] <- path_
  rm(i, path_)
}

# double check file paths
sum(comparison_matrix$sample_a %in% files)
sum(comparison_matrix$sample_b %in% files)

saveRDS(comparison_matrix, paste0("/comparison_matrix_", gsub("-", "_", Sys.Date()), ".rds"))

rm(list = ls())



