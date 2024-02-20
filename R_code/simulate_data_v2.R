library(ggplot2)
library(cowplot)
library(scales)
library(dplyr)
library(tidyr)

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
                           perfectly_diverse = rep(10, 10000),
                           test = c(70000, rep(30000/9999, 9999))[1:10000],
                           test = c(35000, 35000, rep(30000/9999, 9999))[1:10000],
                           test = c(35000, 17500, 17500, rep(30000/9999, 9999))[1:10000],
                           test = c(35000, 11600, 11600, 11600, rep(30000/9999, 9999))[1:10000],
                           test = c(35000, 8750, 8750, 8750, 8750, rep(30000/9999, 9999))[1:10000],
                           test = c(35000, 7000, 7000, 7000, 7000, 7000, rep(30000/9999, 9999))[1:10000],
                           test = c(35000, 5833, 5833, 5833, 5833, 5833, 5833, rep(30000/9999, 9999))[1:10000],
                           test = c(35000, rep(35000/10, 10), rep(30000/9999, 9999))[1:10000],
                           test = c(35000, rep(35000/20, 20), rep(30000/9999, 9999))[1:10000],
                           test = c(35000, rep(35000/50, 50), rep(30000/9999, 9999))[1:10000],
                           test = c(35000, rep(35000/100, 100), rep(30000/9999, 9999))[1:10000]
                           )

sim_tcr_data <- data.frame(clone = paste("clone_", 1:10000),
                           near_perfectly_clonal = c(99999, rep(1/999, 9999)),
                           perfectly_diverse = rep(10, 10000),
                           dominant_clone = c(70000, rep(30000/9999, 9999))[1:10000],
                           oligoclonal_1 = c(35000, 35000, rep(30000/9999, 9999))[1:10000],
                           oligoclonal_2 = c(35000, 17500, 17500, rep(30000/9999, 9999))[1:10000],
                           oligoclonal_3 = c(35000, 11600, 11600, 11600, rep(30000/9999, 9999))[1:10000],
                           oligoclonal_4 = c(35000, 8750, 8750, 8750, 8750, rep(30000/9999, 9999))[1:10000],
                           oligoclonal_5 = c(35000, 7000, 7000, 7000, 7000, 7000, rep(30000/9999, 9999))[1:10000],
                           oligoclonal_6 = c(35000, 5833, 5833, 5833, 5833, 5833, 5833, rep(30000/9999, 9999))[1:10000],
                           oligoclonal_10 = c(35000, rep(35000/10, 10), rep(30000/9999, 9999))[1:10000],
                           oligoclonal_20 = c(35000, rep(35000/20, 20), rep(30000/9999, 9999))[1:10000],
                           oligoclonal_50 = c(35000, rep(35000/50, 50), rep(30000/9999, 9999))[1:10000],
                           oligoclonal_100 = c(35000, rep(35000/100, 100), rep(30000/9999, 9999))[1:10000],
                           oligoclonal_5000 = c(35000, rep(35000/5000, 5000), rep(30000/4999, 4999))[1:10000],
                           oligoclonal_7000 = c(35000, rep(35000/7000, 7000), rep(30000/2999, 2999))[1:10000],
                           oligoclonal_9999 = c(35000, rep(65000/9999, 9999))[1:10000]
)


sim_tcr_data <- data.frame(clone = paste("clone_", 1:10000),
                           near_perfectly_clonal = c(99999, rep(1/999, 9999)),
                           perfectly_diverse = rep(10, 10000),
                           dominant_clone = c(70000, rep(30000/9999, 9999))[1:10000],
                           oligoclonal = c(35000, 35000, rep(30000/9999, 9999))[1:10000],
                           oligoclonal_5 = c(35000, 7000, 7000, 7000, 7000, 7000, rep(30000/9999, 9999))[1:10000],
                           oligoclonal_10 = c(35000, rep(35000/10, 10), rep(30000/9999, 9999))[1:10000],
                           oligoclonal_20 = c(35000, rep(35000/20, 20), rep(30000/9999, 9999))[1:10000],
                           oligoclonal_50 = c(35000, rep(35000/50, 50), rep(30000/9999, 9999))[1:10000],
                           oligoclonal_100 = c(35000, rep(35000/100, 100), rep(30000/9999, 9999))[1:10000],
                           oligoclonal_1000 = c(35000, rep(35000/1000, 1000), rep(30000/8999, 8999))[1:10000],
                           oligoclonal_5000 = c(35000, rep(35000/5000, 5000), rep(30000/4999, 4999))[1:10000],
                           oligoclonal_9999 = c(35000, rep(65000/9999, 9999))[1:10000]
)


## evenness
sample_names <- names(sim_tcr_data)[2:ncol(sim_tcr_data)]
total_clones <- 10000
total_count <- 100000

cumulative_plot <- list()
stats <- data.frame()
for(i in sample_names){
  # clonality
  freq_ <- sim_tcr_data[, i]/100000
  ### Shannon ---
  p_ <- freq_/sum(freq_)
  p_ <- p_[p_ > 0]
  sh_en_ <- -sum(p_ * log2(p_)) # calculate shannon entropy
  sh_cl_ <- 1 - (sh_en_ / log2(length(freq_))) # calculate shannon clonality
  ### Simpson ---
  si_cl_ <- sqrt(sum((freq_ / sum(freq_)) ^ 2))
  ### Gini ---
  freq_g_ <- freq_[order(freq_)]
  cumulative_sum_clones <- seq(1, length(freq_g_))
  perfect_equality_line <- seq(0, 1, length.out = length(freq_g_))
  cumulative_sum_count <- c()
  r_c_s <- 0 # running cumulative sum
  for(j in 1:length(freq_g_)){
    r_c_s <- r_c_s + freq_g_[j]
    cumulative_sum_count <- c(cumulative_sum_count, r_c_s)
    rm(j)
  }
  a_ <- sum(perfect_equality_line * cumulative_sum_clones) - sum(cumulative_sum_count * cumulative_sum_clones)
  g_ <- a_ / sum(perfect_equality_line * cumulative_sum_clones)
  
  stats <- rbind(stats,
                  data.frame(description = i,
                             shannon_entropy = sh_en_,
                             shannon_clonality = sh_cl_,
                             simpson_clonality = si_cl_,
                             gini_coefficient = g_))
  
  # plot
  plot_df <- data.frame(clone = factor(paste("clone_", 1:10000), levels = paste("clone_", 1:10000)), cumulative_sum_clones, freq_)
  to_color_ <- sum(plot_df$freq_ > 0.003)
  if(to_color_ == 0){to_color_ <- 1}

  cumulative_plot[[i]] <- ggplot(plot_df, aes(x = "", y = freq_, fill = clone)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = c(ggplotColours(to_color_), ggplotColours(10000 - to_color_))) +
    coord_polar("y", start=0) +
    theme_void() +
    ggtitle(paste0(i)) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 15))
  
  rm(i, r_c_s)
  rm(list = ls()[grepl("_$", ls())])
}

stats
stats2 <- stats %>%
  pivot_longer(cols = c("shannon_clonality", "simpson_clonality", "gini_coefficient"),
               names_to = "summary_method",
               values_to = "value") %>%
  mutate(description = factor(description,
                              levels = c("near_perfectly_clonal",
                                         "dominant_clone",
                                         "oligoclonal",
                                         paste0("oligoclonal_", c(5, 10, 20, 50, 100, 1000, 5000, 9999)),
                                         "perfectly_diverse")))

ggplot(stats2,
       aes(description, value)) +
  geom_line(aes(group = summary_method),
            linewidth = 2,
            color = "black") +
  geom_line(aes(group = summary_method,
                color = summary_method),
            linewidth = 1.5) +
  scale_color_manual(values = c("hotpink", "dodgerblue", "gold")) +
  guides(colour = guide_legend(ncol = 1, title.position = "top")) +
  theme(axis.text.x = element_text(angle = 60,
                                   hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        legend.position = c(0.25, 0.2),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))

for(i in names(cumulative_plot)){
  print(i)
  jpeg(paste0("~/Documents/BFX_proj/TCRseq_pipeline/_output/simulated_data/", i, ".jpg"), width = 250, height = 275)
  print(cumulative_plot[[i]])
  dev.off()
}


