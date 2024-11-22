# TCRseq tool to mimic most Adaptive analysis while adding additional functionality
# Try as much as possible to write in base R
# try to make universal to all TCRseq inputs
## likely only need nucleic acid sequence, amino acid sequence, count
## using export version 2 from immunoseq Analyzer from Adaptive

# to do
# standard sample summary data [x]
# paired sample comparisons [ ]
# exclusion of highly public clones [ ]

################################################################################
#
# set summary parameters
#
################################################################################

setwd("~/Documents/BFX_proj/TCRseq_pipeline/")

### directories ---
input_dir <- "_input/" # raw data
count_dir <- "_input/counts/" # counts data
output_dir <- "_output/" # target directory for analysis outputs

### identify relevant, TCR-counts-table columns ---
#foo <- read.table(paste0(input_dir, "P9-PB1.tsv"), header = T, sep = "\t")
sep = "\t" # table separation type
nuc_a <- "nucleotide" # nucleic acid sequence
ami_a <- "aminoAcid" # amino acid sequence
count <- "count..templates.reads." # counts columns

### clone filtering parameter ---
# how/whether to condense clone data: "all" = no filtering, "productive" = productive only, "clonotypes" = condense clones to shared AA sequences
clone_group <- "productive" # default is "all"
n_p <- "" # value for non-productive TCR clones in amino acid sequence column
# number of counts to downsample to for the control different read depth
downsample_counts <- 10000

################################################################################
#
# load metadata
#
################################################################################

### ENSURE THERE ARE NO DUPLICATED COLUMN NAMES ###
### ENSURE FILE NAME COLUMN MATCHES ACTUAL FILES AND EXTENSIONS ###

meta <- read.csv(paste0(input_dir, "Chow_PNAS_2020_short.csv")) # load metadata table

################################################################################
#
# load and process individual sample files
#
################################################################################

### initialize data frames for clone and summary data ---
clone_data <- data.frame()
sample_stats <- data.frame()

### initialize data frame for clone data ---
clone_files <- list.files(path = count_dir) # point to target directory

### loop for sample data processing ---
for (file_ in clone_files){
  ##############################################################################
  ### load file ---
  ##############################################################################
  cat(sprintf("Loading '%s'\n", file_)) # print progress
  data_ <- read.table(paste0(count_dir, file_), header = T, sep = sep) # load data
  #$ possible error where file is missing key column names
  data_ <- data.frame(row.names = data_[, nuc_a],
                      nucleic_acid = data_[, nuc_a],
                      amino_acid = data_[, ami_a],
                      count = data_[, count],
                      file = file_)
  
  ##############################################################################
  ### filter clones ---
  ##############################################################################
  if (clone_group == "all"){
    data_ <- data_
  } else if (clone_group == "productive") {
    data_ <- data_[data_$amino_acid != n_p, ] # filter out out-of-frame clones
  } else if (clone_group == "clonotypes"){
    data_ <- data_[data_$amino_acid != n_p, ] # filter out out-of-frame clones
    data_ <- data_[order(data_$count, decreasing = T), ] # put dominant nucleic acid sequences on top
    data__ <- data_[F, ] # initialize new data table to populate with clonotype reduced data
    ### for loop to condense clonotype data
    for(a_a_ in unique(data_$amino_acid)){
      data___ <- data_[data_$amino_acid == a_a_, ]
      data___ <- data.frame(nucleic_acid = data___$nucleic_acid[1], # dominant nucleic acid sequence
                            amino_acid = data___$amino_acid[1], # shared amino acid sequence
                            count = sum(data___$count), # sum of all counts
                            file = data___$file[1]) # file name
      data__ <- rbind(data__, data___)
      rm(a_a_, data___) # clean up
    }
    data_ <- data__ # overwrite uncondensed data frame
    rm(data__) # clean up
  } else {
    data_ <- data_
  }
  
  ##############################################################################
  ### downsample data ---
  ##############################################################################
  if (is.numeric(downsample_counts)){
    ### downsampling ---
    cat(sprintf("Sample: %s, subsampling %d to %d random counts\n", file_, sum(data_$count), downsample_counts))
    data_old_ <- cbind(data_[, c("nucleic_acid", "count")]) # extract parent counts data
    data_old_ <- rep(data_old_$nucleic_acid, times = data_old_$count) # convert to vector, repeat each template for number of counts
    set.seed(415); data_new_ <- data_old_[sample.int(length(data_old_), length(data_old_))] # shuffle counts (probably doesn't matter)
    set.seed(415); data_new_ <- data_new_[sample.int(length(data_new_), downsample_counts)] # subsample counts
    data_new_ <- data.frame(table(nucleic_acid = data_new_)); rownames(data_new_) <- data_new_$nucleic_acid # convert to data.frame
    
    ### amend counts ---
    data_ <- data_[data_$nucleic_acid %in% data_new_$nucleic_acid, ] # remove data subsampled out
    data_new_ <- data_new_[data_$nucleic_acid, ] # rearrange to match data_
    data_$count <- data_new_$Freq
    
    rm(data_new_) # clean up
  }
  
  ##############################################################################
  ### calculate frequency ---
  ##############################################################################
  data_$freq <- data_$count/sum(data_$count)
  data_ <- data_[, c("file", "nucleic_acid", "amino_acid", "count", "freq")]
  
  ##############################################################################
  ### add sample data to clone data ---
  ##############################################################################
  clone_data <- rbind(clone_data, data_)
  
  ##############################################################################
  ### summarize sample stats ---
  ##############################################################################
  ### gross features ---
  stats_ <- data.frame(
    row.names = file_,
    file = file_,
    unique_clones = nrow(data_),
    total_count = sum(data_$count),
    count_min = min(data_$count),
    count_max = max(data_$count),
    count_mean = mean(data_$count),
    freq_min = min(data_$freq),
    freq_max = max(data_$freq),
    freq_mean = mean(data_$freq)
    )
  
  ### clonality ---
  freq_ <- data_$freq
  ### Shannon ---
  p_ <- freq_/sum(freq_)
  p_ <- p_[p_ > 0]
  sh_en_ <- -sum(p_ * log2(p_)) # calculate shannon entropy
  sh_cl_ <- 1 - (sh_en_ / log2(length(freq_))) # calculate shannon clonality
  ### Simpson ---
  si_cl_ <- sqrt(sum((freq_ / sum(freq_)) ^ 2))
  ### append sample stats ---
  stats_ <- cbind(stats_,
                  data.frame(shannon_entropy = sh_en_,
                             shannon_clonality = sh_cl_,
                             simpson_clonality = si_cl_))
  
  rm(freq_, p_, sh_en_, sh_cl_, si_cl_) # clean up
  
  ### clone size groups ---
  cl_gp_ <- c("XL", # (>1%)",
              "L", # (1-0.1%)",
              "M", # (0.1-0.01%)", 
              "S", # (0.01-0.001%)",
              "XS") # (<0.001%)")
  
  data_$clone_group <- ifelse(data_$freq > 0.01, cl_gp_[1],
                              ifelse(data_$freq > 0.001, cl_gp_[2],
                                     ifelse(data_$freq > 1E-04, cl_gp_[3],
                                            ifelse(data_$freq > 1E-05, cl_gp_[4],
                                                   ifelse(data_$freq <= 1E-05, cl_gp_[5])))))
  data_$clone_group <- factor(data_$clone_group, levels = cl_gp_)
  
  cg_stats_ <- data.frame(matrix(ncol = 15, nrow = 1,
                                 dimnames = list(file_, c(paste0(cl_gp_[1:5], "_clone_no"),
                                                          paste0(cl_gp_[1:5], "_clone_freq"),
                                                          paste0(cl_gp_[1:5], "_rep_freq")))))
  for(cg_ in cl_gp_){
    data_cg_ <- data_[data_$clone_group == cg_, ]
    cg_stats_[file_, paste0(cg_, "_clone_no")] <- nrow(data_cg_)
    cg_stats_[file_, paste0(cg_, "_clone_freq")] <- signif(nrow(data_cg_)/nrow(data_), 6)
    cg_stats_[file_, paste0(cg_, "_rep_freq")] <- signif(sum(data_cg_$count)/sum(data_$count), 6)
    rm(cg_) # clean up
  }
  ### append sample stats ---
  stats_ <- cbind(stats_,
                  cg_stats_)
  
  ##############################################################################
  ### add sample data to sample summary data ---
  ##############################################################################
  sample_stats <- rbind(sample_stats, stats_)
  
  rm(list = ls()[grepl("_$", ls())]) # clean up
}

write.csv(clone_data, paste0(output_dir, "clone_data.csv"))
write.csv(sample_stats, paste0(output_dir, "sample_summary.csv"))




