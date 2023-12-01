# Additional diversity parameter calculator
## run D50 and Di index calculation

additional_diversity <- data.frame()

for(i in unique(clone_data$file)[1]){
  data_ <- clone_data[clone_data$file == i, ] # subset 
  data_ <- data_[order(data_$count, decreasing = T), ] # order clones by count
  
  # take top 10000 by rank sum
  if(nrow(data_) > 10000){
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
  
  data_ <- data.frame(data_,
                      cumulative_sum_clones = cumulative_sum_clones,
                      cumulative_sum_count = cumulative_sum_count/total_count*100)
  
  auc <- data_
  
}
