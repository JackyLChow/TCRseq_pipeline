# TCRseq tool to mimic most Adaptive analysis while adding additional functionality
# Try as much as possible to write in base R
# try to make universal to all TCRseq inputs
## likely only need nucleic acid sequence, amino acid sequence, count
## using export version 2 from immunoseq Analyzer from Adaptive

# to do
# standard sample summary data [x]
# paired sample comparisons [x]
# exclusion of highly public clones [ ]

### remember to run separate parameters script

# v_large:
## When risk of clone_data being too large, write each condensed file to RDS for compilation/referencing later
# v_2
## Break pipeline into two parts:
### 1) Data munching: condense clone data, filter top 10k clones AND filter top 10k productive clones
### 2) Data summarizing

################################################################################
#
# load and process individual sample files
#
################################################################################
load("")

processed <- 0

clone_files <- files[1:1000]

start_time <- Sys.time()
### loop for sample data processing ---
for (file_ in clone_files){
  ##############################################################################
  ### load file ---
  ##############################################################################
  file_time_ <- Sys.time()
  data_ <- read.table(paste0(file_), header = T, sep = sep) # load data
  # trim file path
  file_ <- unlist(strsplit(file_, "/"))
  file_ <- file_[grepl("airr.csv$", file_)]
  cat(sprintf("Loading '%s'\n", file_)) # print progress
  
  if(!is.null(keep_cols)){
    keep_data_ <- data_[, keep_cols]
  }
  data_ <- data.frame(#row.names = data_[, nuc_a],
    nucleic_acid = data_[, nuc_a],
    amino_acid = data_[, ami_a],
    count = data_[, count],
    file = file_,
    clone_id = data_[, cln_id])
  if(!is.null(keep_cols)){
    data_ <- data.frame(data_, keep_data_)
  }
  
  ##############################################################################
  ### condense clones ---
  ##############################################################################
  if (clone_group == "all"){
    cat("No clone condensing\n") # print progress
    data_ <- data_
  } else if (clone_group %in% colnames(data_)){ # custom column for clone id
    cat(sprintf("Condensing clones by '%s'\n", clone_group)) # print progress
    data_ <- data_[order(data_$count, decreasing = T), ] # put dominant nucleic acid sequences on top
    data__ <- data_[F, ] # initialize new data table to populate with clonotype reduced data
    ### for loop to condense clonotype data
    for(id_ in unique(data_[, clone_group])){
      data___ <- data_[data_[, clone_group] == id_, ]
      if(!is.null(keep_cols)){
        keep_data___ <- data___[1, keep_cols]
      }
      data___ <- data.frame(nucleic_acid = data___$nucleic_acid[1], # dominant nucleic acid sequence
                            amino_acid = data___$amino_acid[1], # shared amino acid sequence
                            count = sum(data___$count), # sum of all counts
                            file = data___$file[1],
                            clone_id = data___$clone_id[1]) # file name
      if(!is.null(keep_cols)){
        data___ <- data.frame(data___, keep_data___)
      }
      data__ <- rbind(data__, data___)
      rm(data___) # clean up
    }
    data_ <- data__ # overwrite un-condensed data frame
    rm(data__) # clean up
  } else {
    cat("No clone condensing\n") # print progress
    data_ <- data_
  }
  
  
  
  ##############################################################################
  ### downsample data ---
  ##############################################################################
  # all clones
  if (is.numeric(rank_cutoff)){
    if (nrow(data_) > rank_cutoff){
      ### trim by rank dominance ---
      cat(sprintf("Sample: %s, downsampling %d to top %d clones\n", file_, nrow(data_), rank_cutoff))
      data__ <- data_[order(data_$count, decreasing = T), ] # order by dominance
      data__ <- data__[1:rank_cutoff, ] # apply cutoff
    }
    if (nrow(data_) < rank_cutoff){
      cat(sprintf("Sample: %s, cannot downsample %d to top %d clones\n", file_, nrow(data_), rank_cutoff))
      data__ <- data_
    }
    data__$freq <- data__$count/sum(data__$count)
    
    saveRDS(data__[, c("file", keep_cols, "nucleic_acid", "amino_acid", "clone_id", "count", "freq")], paste0(output_dir, "/top_10k_all/", gsub("\\.", "_", file_), ".rds"))
    # glimpse(data__)
    rm(data__) # clean up
  }
  
  # productive only
  cat(paste("Removing", sum(grepl(s_c, data_$amino_acid) | data_$amino_acid == n_p), "clones\n"))
  data_ <- data_[!grepl(s_c, data_$amino_acid) & data_$amino_acid != n_p, ] # filter out unproductive clones
  if (is.numeric(rank_cutoff)){
    if (nrow(data_) > rank_cutoff){
      ### trim by rank dominance ---
      cat(sprintf("Sample: %s, downsampling %d to top %d clones\n", file_, nrow(data_), rank_cutoff))
      data__ <- data_[order(data_$count, decreasing = T), ] # order by dominance
      data__ <- data__[1:rank_cutoff, ] # apply cutoff
    }
    if (nrow(data_) < rank_cutoff){
      cat(sprintf("Sample: %s, cannot downsample %d to top %d clones\n", file_, nrow(data_), rank_cutoff))
      data__ <- data_
    }
    data__$freq <- data__$count/sum(data__$count)
    
    saveRDS(data__[, c("file", keep_cols, "nucleic_acid", "amino_acid", "clone_id", "count", "freq")], paste0(output_dir, "/top_10k_productive/", gsub("\\.", "_", file_), ".rds"))
    # glimpse(data__)
  }
  
  processed <- processed + 1
  
  cat(paste("This file processing time:\n"))
  print(Sys.time() - file_time_)
  
  cat(paste("Total samples from this cohort processed:", processed, "\n"))
  
  cat(paste("Total samples processed:",
            length(list.files(output_dir, recursive = T))/2,
            "\n"))
  print(Sys.time() - start_time)
  cat("\n")
}
