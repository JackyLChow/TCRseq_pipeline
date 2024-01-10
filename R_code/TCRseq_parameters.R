# TCRseq parameters
library(dplyr)

adaptive_out <- read.table("~/Documents/BFX_proj/TCRseq_pipeline/_input/Chow_PNAS_2020_Adaptive_summary.tsv", header = T, sep = "\t")

### directories ---
count_dir <- "~/Documents/BFX_proj/TCRseq_pipeline/_input/counts/" # counts data
output_dir <- "~/Documents/BFX_proj/TCRseq_pipeline/_output/" # target directory for analysis outputs
output_dir_temp <- "~/Documents/BFX_proj/TCRseq_pipeline/_output/temp/" # target directory for analysis outputs

### identify relevant, TCR-counts-table columns ---
sep = "\t" # table separation type
nuc_a <- "nucleotide" # nucleic acid sequence
ami_a <- "aminoAcid" # amino acid sequence
count <- "count..templates.reads." # counts columns
cln_id <- "nucleotide"
keep_cols <- c("vGeneName", "jGeneName")

### clone filtering parameter ---
# how/whether to condense clone data: "all" = no filtering, "productive" = productive only, "clonotypes" = condense clones to shared AA sequences
clone_group <- "all" # default is "all"; "productive" is only productive sequences"; "aa" is only productive sequences, condensed by common aa sequence
n_p <- "" # value for non-productive TCR clones in amino acid sequence column
s_c <- "\\*" # value to grep for early stop codons
resample_counts <- NULL # number of counts to resample to for the control different read depth
rank_cutoff <- 10000L # number of counts to resample to for the control different read depth