library(ggplot2)

# Generate simulated TCRseq data to explore different diversity metrics; 100k total counts, 10k clones

## simulated dataframe
sim_tcr_data <- data.frame(clone = paste("clone_", 1:10000),
                           near_perfectly_clonal = c(99999, rep(1/999, 9999)),
                           highly_clonal = c(40000, rep(60000/9999, 9999)),
                           oligoclonal = c(40000, 20000, 10000, rep(30000/9997, 9997)),
                           perfectly_diverse = rep(10, 10000))

## evenness
sample_names <- names(sim_tcr_data)[2:5]
total_clones <- 10000
total_count <- 100000

cumulative_plot <- list()

for(i in sample_names){
  # clonality
  freq_ <- sim_tcr_data[, i]/100000
  ## Shannon
  p_ <- freq_/sum(freq_)
  p_ <- p_[p_ > 0]
  sh_en_ <- -sum(p_ * log2(p_)) # calculate shannon entropy
  sh_cl_ <- 1 - (sh_en_ / log2(length(freq_))) # calculate shannon clonality
  ## Simpson
  si_cl_ <- sqrt(sum((freq_ / sum(freq_)) ^ 2))
  
  ## cumulative data
  cumulative_sum_clones <- seq(1, 100, length.out = total_clones)
  cumulative_sum_freq <- c()
  r_c_s <- 0 # running cumulative sum
  for(j in 1:nrow(sim_tcr_data)){
    r_c_s <- r_c_s + (sim_tcr_data[j, i] / total_count * 100)
    cumulative_sum_freq <- c(cumulative_sum_freq, r_c_s)
  }
  
  d50 <- match(T, cumulative_sum_freq >= 50)/100
  d_i <- (10000 - sum(cumulative_sum_freq * (cumulative_sum_clones[2] - cumulative_sum_clones[1])))/100
  
  # plot
  plot_df <- data.frame(id = as.character(1:10000), cumulative_sum_clones, cumulative_sum_freq)
  
  ggplot(data, aes(x="", y=value, fill=group)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0)
  
  ggplot(plot_df, aes(x = "lol", y = cumulative_sum_freq, fill = id)) +
    geom_bar(position = "stack", stat = "identity")
  ggplot(plot_df, aes(cumulative_sum_clones, cumulative_sum_freq)) +
    geom_area(alpha = 0.5) +
    geom_vline(xintercept = 50, color = "red", linetype = "dotted") +
    geom_segment(x = 0, y = 50, xend = d50, yend = 50, color = "red") +
    geom_segment(x = d50, y = 50, xend = d50, yend = 0, color = "red") +
    theme(aspect.ratio = 1)
  
  print(paste(i, sh_en_, sh_cl_, si_cl_, d50, d_i))
}

