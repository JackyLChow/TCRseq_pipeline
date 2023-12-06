library(ggplot2)
library(cowplot)
library(scales)

# Generate simulated TCRseq data to explore different diversity metrics; 100k total counts, 10k clones

## ggplot customization
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

## simulated dataframe
sim_tcr_data <- data.frame(clone = paste("clone_", 1:10000),
                           near_perfectly_clonal = c(99999, rep(1/999, 9999)),
                           dominant_clone = c(70000, rep(30000/9999, 9999)),
                           oligoclonal_1 = c(40000, 20000, 10000, rep(30000/9997, 9997)),
                           oligoclonal_2 = c(40000, 10000, 10000, 5000, 5000, rep(30000/9993, 9995)),
                           singularly_clonal = c(40000, rep(60000/9999, 9999)),
                           perfectly_diverse = rep(10, 10000))

## evenness
sample_names <- names(sim_tcr_data)[2:7]
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
  plot_df <- data.frame(clone = factor(paste("clone_", 1:10000), levels = paste("clone_", 1:10000)), cumulative_sum_clones, freq_)
  
  cumulative_plot[[i]] <- plot_grid(
    ggplot(plot_df, aes(x = "", y = freq_, fill = clone)) +
      geom_bar(stat = "identity", width = 1) +
      scale_fill_manual(values = c(ggplotColours(6), ggplotColours(9994))) +
      coord_polar("y", start=0) +
      theme_void() +
      ggtitle(paste0(i)) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5))
    ,
    ggplot(plot_df, aes(cumulative_sum_clones, cumulative_sum_freq)) +
      ggtitle(paste0("Shannon entropy: ", round(sh_en_, 3),
                     "\nShannon clonality: ", round(sh_cl_, 3),
                     "\nSimpson clonality: ", round(si_cl_, 3),
                     "\nD50: ", round(d50, 3),
                     "\ndiversity index: ", round(d_i, 1))) +
      xlab("% of richness") +
      ylab("% of repertoire") +
      geom_area(alpha = 0.5) +
      geom_vline(xintercept = 50, color = "red", linetype = "dotted") +
      geom_segment(x = 0, y = 50, xend = d50, yend = 50, color = "red") +
      geom_segment(x = d50, y = 50, xend = d50, yend = 0, color = "red") +
      theme(aspect.ratio = 1) 
    , nrow = 2
  )
  print(paste(i, round(sh_en_, 3), round(sh_cl_, 3), round(si_cl_, 3), round(d50, 3), round(d_i, 1)))
}

plot_grid(plotlist = cumulative_plot, nrow = 1)

