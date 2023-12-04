# TCRseq_differential_abundance
## For each clone, compare relative abundances

clone_data <- readRDS("~/Documents/BFX_proj/TCRseq_pipeline/_output/clone_data.rds")
differential_clone_abundance_results_folder <- "~/Documents/BFX_proj/TCRseq_pipeline/_output/differential_clone_abundance/"

comparison_matrix <- data.frame(sample_a = "P10-PB1.tsv",
                                sample_b = unique(clone_data$file),
                                id = "nucleic_acid",
                                d_a_method = "fisher")

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
    
    counts_a <- clone_dat[clone_dat$file == a_, c(id_, "count")][1:10, ]
    if(any(duplicated(counts_a[, id_]))){
      stop("duplicated ids in counts_A")
    }
    counts_b <- clone_dat[clone_dat$file == b_, c(id_, "count")][1:10, ]
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

# differential clone abundance summarizer
dca_files <- list.files(differential_clone_abundance_results_folder, full.names = T)

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
                       contracted = sum(res_$significance == "A > A")
    )
    dca_summary <- rbind(dca_summary, res_)
  }
  rm(i)
  return(dca_summary)
}

foo <- differential_clone_abundance_summarizer()










