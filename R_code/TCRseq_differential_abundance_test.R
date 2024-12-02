# TCRseq_differential_abundance
## For each clone, compare relative abundances

###
### binomial and fisher exact need integers, must use counts or umi
### fisher and binomial sensitive to read depth
###

differential_clone_abundance_calculator <- function(
    clone_dat = clone_data,
    comparison_mtx = comparison_matrix
){
  # prepare output data
  comparison_results <- data.frame()
  # run comparisons
  for(i in 1:nrow(comparison_mtx)){
    cat(paste("running", comparison_mtx[i, "sample_a"], "vs", comparison_mtx[i, "sample_b"], "by", comparison_mtx[i, "id"], "\n"))
    
    a_ <- comparison_mtx[i, "sample_a"]
    b_ <- comparison_mtx[i, "sample_b"]
    id_ <- comparison_mtx[i, "id"]
    d_a_m_ <- comparison_mtx[i, "d_a_method"]
    
    counts_a <- clone_dat[clone_dat$file == a_, c(id_, "count")]
    counts_a$count <- counts_a$count
    if(any(duplicated(counts_a[, id_]))){
      stop("duplicated ids in counts_A")
    }
    counts_b <- clone_dat[clone_dat$file == b_, c(id_, "count")]
    counts_a$count <- counts_a$count
    if(any(duplicated(counts_b[, id_]))){
      stop("duplicated ids in counts_B")
    }
    
    counts <- merge(counts_a, counts_b, by = id_, all = TRUE, suffixes = c("_a", "_b"))
    counts[is.na(counts)] <- 0
    
    A <- counts$count_a
    B <- counts$count_b
    sumA <- sum(A)
    sumB <- sum(B)
    
    # Differential clone abundance
    p_value <- c()
    for(j in 1:nrow(counts)){
      if(d_a_m_ == "binomial"){
        p_value <- c(p_value,
                     binom.test(counts[j, "count_a"], counts[j, "count_a"] + counts[j, "count_b"], p = sumA/(sumA + sumB), alternative = "two.sided")$p.value) # p_value for binomial test
      }
      if(d_a_m_ == "fisher"){
        p_value <- c(p_value,
                     fisher.test(matrix(c(counts[j, "count_a"], sumA, counts[j, "count_b"], sumB), nrow = 2), alternative = "two.sided")$p.value)
      }
    }
    
    differential_abundance <- data.frame(sample_a = a_, sample_b = b_,
                                         id = counts[, 1], counts[, c("count_a", "count_b")],
                                         p_value = p_value,
                                         p_adj = p.adjust(p_value, method = "BH"),
                                         p_value_method = d_a_m_)
    differential_abundance$freq_a <- differential_abundance$count_a/sum(differential_abundance$count_a)
    differential_abundance$freq_b <- differential_abundance$count_b/sum(differential_abundance$count_b)
    differential_abundance$significance <- ifelse(
      differential_abundance$p_adj <= 0.01 & differential_abundance$freq_a > differential_abundance$freq_b, "A > B",
      ifelse(
        differential_abundance$p_adj <= 0.01 & differential_abundance$freq_a < differential_abundance$freq_b, "B > A",
        ifelse(
          differential_abundance$p_adj > 0.01, "not_significant", "error"
        )
      )
    )
    saveRDS(differential_abundance, paste0(differential_clone_abundance_results_folder, gsub("\\.", "_", a_), "_vs_", gsub("\\.", "_", b_), ".rds"))
    cat(paste("completed\n"))
  }
  rm(i, j)
}

differential_clone_abundance_calculator()

foo <- readRDS("~/Documents/BFX_proj/TCRseq_test/_output/differential_clone_abundance/801_4C_New_repeat_2308_RD_0021_YPL5_D501_D701_SER1_BC801_TRA_2308_RD_0021_YPL5_D501_D701_SER1_BC801_TRA_pep_vs_801_4C_New_2308_RD_0021_YPL5_A501_A701_SER1_BC801_TRA_2308_RD_0021_YPL5_A501_A701_SER1_BC801_TRA_pep.rds")
ggplot(foo, aes(log10(freq_a + 1), log10(freq_b + 1))) +
  geom_point()

bar <- readRDS("~/Documents/BFX_proj/TCRseq_test/_output/differential_clone_abundance/801_4C_New_repeat_2308_RD_0021_YPL5_D501_D701_SER1_BC801_TRA_2308_RD_0021_YPL5_D501_D701_SER1_BC801_TRA_pep_vs_811_4C_New_repeat_2308_RD_0021_YPL5_D502_D702_SER1_BC811_TRA_2308_RD_0021_YPL5_D502_D702_SER1_BC811_TRA_pep.rds")
ggplot(bar, aes(log10(freq_a + 1), log10(freq_b + 1))) +
  geom_point()

# differential clone abundance summarizer
dca_files <- list.files(differential_clone_abundance_results_folder, full.names = T, pattern = ".rds")

differential_clone_abundance_summarizer <- function(
    results_files = dca_files){
  dca_summary <- data.frame()
  for(i in results_files){
    res_ <- readRDS(i)
    res_ <- data.frame(sample_a = unique(res_$sample_a),
                       sample_b = unique(res_$sample_b),
                       total = nrow(res_),
                       significant = sum(res_$significance != "not_significant"),
                       not_significant = sum(res_$significance == "not_significant"),
                       expanded = sum(res_$significance == "B > A"),
                       contracted = sum(res_$significance == "A > B")
    )
    dca_summary <- rbind(dca_summary, res_)
  }
  rm(i)
  return(dca_summary)
}

differential_clone_abundance_summary <- differential_clone_abundance_summarizer()

rm(dca_files)




