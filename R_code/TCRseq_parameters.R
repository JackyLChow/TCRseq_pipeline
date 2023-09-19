# TCRseq parameters

### directories ---
count_dir <- "~/Documents/BFX_proj/TCRseq_pipeline/_input/counts/" # counts data
output_dir <- "~/Documents/BFX_proj/TCRseq_pipeline/_output/" # target directory for analysis outputs

### identify relevant, TCR-counts-table columns ---
sep = "\t" # table separation type
nuc_a <- "nucleotide" # nucleic acid sequence
ami_a <- "aminoAcid" # amino acid sequence
count <- "count..templates.reads." # counts columns

### clone filtering parameter ---
# how/whether to condense clone data: "all" = no filtering, "productive" = productive only, "clonotypes" = condense clones to shared AA sequences
clone_group <- "all" # default is "all"
n_p <- "" # value for non-productive TCR clones in amino acid sequence column
# number of counts to downsample to for the control different read depth
downsample_counts <- NULL