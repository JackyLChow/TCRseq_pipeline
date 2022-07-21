# TCRseq tool to mimic most Adaptive analysis while adding additional functionality
# Try as much as possible to write in base R
# try to make universal to all TCRseq inputs
## likely only need nucleic acid sequence, amino acid sequence, count
## using export version 2 from immunoseq Analyzer from Adaptive

# to do
# standard sample summary data [ ]
# paired sample comparisons [ ]
# exclusion of highly public clones [ ]

#################################
## Use Adaptive v1 export only ##
#################################

## Load packages ##
# CRAN
#library(dplyr)
#library(tidyr)
#library(ggplot2)
#library(readxl)
#library(xlsx)

################################################################################
#
# set summary parameters
#
################################################################################

setwd("~/Documents/BFX_proj/TCRseq_pipeline/")

### directories ---
raw_dir <- "data/raw/" # raw data
count_dir <- "data/raw/counts/" # counts data
der_dir <- "data/derived/" # target directory for analysis outputs

### identify relevant, TCR-counts-table columns ---
#foo <- read.table(paste0(raw_dir, "P9-PB1.tsv"), header = T, sep = "\t")
sep = "\t" # table separation type
nuc_a <- "nucleotide" # nucleic acid sequence
ami_a <- "aminoAcid" # amino acid sequence
count <- "count..templates.reads." # counts columns

### clone filtering parameter ---
# how/whether to condense clone data: "all" = no filtering, "productive" = productive only, "clonotypes" = condense clones to shared AA sequences
clone_group <- "productive" # default is "all"
# number of counts to subsample to for the control different read depth
subsample <- 10000

################################################################################
#
# load metadata
#
################################################################################

### ENSURE THERE ARE NO DUPLICATED COLUMN NAMES ###
### ENSURE FILE NAME COLUMN MATCHES ACTUAL FILES AND EXTENSIONS ###

meta <- read.csv(paste0(raw_dir, "Chow_PNAS_2020_short.csv")) # load metadata table

################################################################################
#
# load and process individual sample files
#
################################################################################

### initialize data frames for clone and summary data ---
clone_data <- data.frame()
sample_data <- data.frame()

### initialize data frame for clone data ---
clone_files <- list.files(path = count_dir) # point to target directory

### loop for sample data processing ---
for (file_ in clone_files){
  ### load file ---
  cat(sprintf("Loading '%s'\n", file_)) # print progress
  data_ <- read.table(paste0(count_dir, file_), header = T, sep = sep) # load data
  #$ possible error where file is missing key column names
  data_ <- data.frame(row.names = data_[, nuc_a],
                      nucleic_acid = data_[, nuc_a],
                      amino_acid = data_[, ami_a],
                      count = data_[, count],
                      file = file_)
  
  ### filter clones ---
  if (clone_group == "all"){
    data_ <- data_
  } else if (clone_group == "productive") {
    data_ <- data_[data_$amino_acid != "", ] # filter out out-of-frame clones
  } else if (clone_group == "clonotypes"){
    data_ <- data_[data_$amino_acid != "", ] # filter out out-of-frame clones
    data_ <- data_[order(data_$count, decreasing = T), ] # put dominant nucleic acid sequences on top
    data__ <- data.frame() # initialize new data table to populate with clonotype reduced data
    ### for loop to condense clonotype data
    for(a_a_ in unique(data_$amino_acid)){
      data___ <- data_[data_$amino_acid == a_a_, ]
      data___ <- c(nucleic_acid = data___$nucleic_acid[1],
                   amino_acid = data___$amino_acid[1],
                   count = sum(data___$count),
                   file = data___$file[1])
      data__ <- rbind(data__, data___)
      rm(a_a_, data___) # clean up
    }
    data_ <- data__ # overwrite uncondensed data frame
    rm(data__) # clean up
  } else {
    data_ <- data_
  }
  
  ### downsample data ---
  
  
  
  clone_data <- rbind(clone_data, data_)
  
  ### summarize data ---
  
}
