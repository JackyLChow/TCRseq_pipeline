library(DescTools)

# Additional diversity parameter calculator

clone_data <- readRDS("~/Documents/BFX_proj/INT/processed_data/TCRseq/iRep_2023_12/munch_1_clone_data.rds")

additional_diversity <- data.frame()

for(i in unique(clone_data$file)){
  data_ <- clone_data[clone_data$file == i, ] # subset
  data_ <- data_[order(data_$count, decreasing = T), ] # order clones by count
  
  # take top 10000 by rank sum
  if(nrow(data_) > 10000){
    cat(paste("\ndownsample to 10000:", i, "\n"))
    data_ <- data_[1:10000, ]
  }

  # sample data paramters
  total_clones <- nrow(data_)
  total_count <- sum(data_$count)
  
  # cumulative data
  cumulative_sum_clones <- seq(1, 100, length.out = total_clones)
  cumulative_sum_count <- c()
  r_c_s <- 0 # running cumulative sum
  for(j in 1:nrow(data_)){
    r_c_s <- r_c_s + data_[j, "count"]
    cumulative_sum_count <- c(cumulative_sum_count, r_c_s)
  }
  cumulative_sum_freq <- cumulative_sum_count/total_count * 100
  
  d50 <- match(T, cumulative_sum_freq >= 50)/total_clones * 100
  diversity_index <- (10000 - sum(cumulative_sum_freq * (cumulative_sum_clones[2] - cumulative_sum_clones[1])))/100
  # diversity_index <- (10000 - AUC(cumulative_sum_clones, cumulative_sum_freq))/100 # using AUC function from DescTools
  
  additional_diversity <- rbind(additional_diversity,
                                data.frame(file = i,
                                           d50 = d50,
                                           diversity_index = diversity_index))
  rm(i, j, cumulative_sum_clones, cumulative_sum_count, cumulative_sum_freq, d50, r_c_s, total_clones, total_count, diversity_index, data_)
}

sample_stats <- left_join(sample_stats, additional_diversity, by = "file")
saveRDS(sample_stats, paste0(output_dir, "sample_stats.rds"))

rm(additional_diversity)
