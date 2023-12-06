# Generate simulated TCRseq data to explore different diversity metrics; 100k total counts, 10k clones

## simulated dataframe
sim_tcr_data <- data.frame(clone = paste("clone_", 1:10000),
                           perfectly_clonal = c(99999, rep(1/999, 9999)),
                           highly_clonal = c(40000, rep(60000/9999, 9999)),
                           oligoclonal = c(40000, 20000, 10000, rep(30000/9997, 9997)),
                           perfectly_diverse = rep(10, 10000))

## evenness
sample_names <- names(sim_tcr_data)[2:5]

for(i in sample_names){
  ### Shannon ---
  freq_ <- sim_tcr_data[, i]/10000
  p_ <- freq_/100000
  p_ <- p_[p_ > 0]
  sh_en_ <- -sum(p_ * log2(p_)) # calculate shannon entropy
  sh_cl_ <- 1 - (sh_en_ / log2(10000)) # calculate shannon clonality
  ### Simpson ---
  si_cl_ <- sqrt(sum((freq_ / 100000) ^ 2))
}

