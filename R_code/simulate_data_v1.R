# Generate simulated TCRseq data to explore different diversity metrics; 100k total counts, 10k clones

## simulated dataframe
sim_tcr_data <- data.frame(clone = paste("clone_", 1:10000),
                           near_perfectly_clonal = c(90001, rep(1, 9999)),
                           highly_clonal = c(40000, rep(60000/9999, 9999)),
                           oligoclonal = c(40000, 20000, 10000, rep(30000/9997, 9997)),
                           perfectly_diverse = rep(10, 10000))

### Shannon ---
p_ <- freq_/sum(freq_)
p_ <- p_[p_ > 0]
sh_en_ <- -sum(p_ * log2(p_)) # calculate shannon entropy
sh_cl_ <- 1 - (sh_en_ / log2(length(freq_))) # calculate shannon clonality
### Simpson ---
si_cl_ <- sqrt(sum((freq_ / sum(freq_)) ^ 2))