# TCRseq_differential_abundance
## For each clone, compare relative abundances

# TCRseq pipeline parameters to match Adaptive: clone is nucleic acid, no filtering, no downsampling
# foo <- read.table("~/Documents/BFX_proj/TCRseq_pipeline/_input/Results_differential-abundance_2023-09-19_03-55/P10-PB1_VS_P10-PB2.differentialAbundance.tsv", sep = "\t", header = T)

# v_large
## pull individual condensed and downsampled files to run comparisons
differential_clone_abundance_results_folder <- 
comparison_matrix <- readRDS()

differential_clone_abundance_calculator <- function(
    comparison_mtx = comparison_matrix){
  start_time <- Sys.time()
  # prepare output data
  comparison_results <- data.frame()
  # run comparisons
  for(i in 1:nrow(comparison_mtx)){
    file_time <- Sys.time()
    cat(paste("running", comparison_mtx[i, "sample_a"], "vs", comparison_mtx[i, "sample_b"], "by", comparison_mtx[i, "id"], "\n"))
    
    a_ <- comparison_mtx[i, "sample_a"]
    b_ <- comparison_mtx[i, "sample_b"]
    bc_a_ <- comparison_mtx[i, "name_a"]
    bc_b_ <- comparison_mtx[i, "name_b"]
    id_ <- comparison_mtx[i, "id"]
    d_a_m_ <- comparison_mtx[i, "d_a_method"]
    
    counts_a <- readRDS(a_)
    counts_a[, "cdr3_v_j"] <- gsub(" ", "_", paste(counts_a$cdr3_aa, counts_a$v_call, counts_a$j_call)) # fix in pipeline
    
    counts_a <- counts_a[, c(id_, "count")]
    if(any(duplicated(counts_a[, id_]))){
      stop("duplicated ids in counts_A")
    }
    counts_b <- readRDS(b_)
    counts_b[, "cdr3_v_j"] <- gsub(" ", "_", paste(counts_b$cdr3_aa, counts_b$v_call, counts_b$j_call)) # fix in pipeline
    
    counts_b <- counts_b[, c(id_, "count")]
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
                                         barcode_chain_a = bc_a_, 
                                         barcode_chain_b = bc_b_,
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
    glimpse(differential_abundance[order(differential_abundance$p_adj), ])
    
    cat(paste("This comparison processing time:\n"))
    print(Sys.time() - file_time)
    
    cat(paste("Total comparisons processed:", i, "\n"))
    print(Sys.time() - start_time)
    
    saveRDS(differential_abundance, paste0(differential_clone_abundance_results_folder, "/", gsub("\\.", "_", bc_a_), "_vs_", gsub("\\.", "_", bc_b_), ".rds"))
    cat(paste("completed\n\n"))
  }
  rm(i, j)
}

differential_clone_abundance_calculator(comparison_mtx = comparison_matrix)

# differential clone abundance summarizer
dca_files <- list.files(differential_clone_abundance_results_folder, full.names = T, pattern = ".rds$")

differential_clone_abundance_summarizer <- function(
    results_files = dca_files){
  dca_summary <- data.frame()
  for(i in results_files){
    res_ <- readRDS(i)
    res_ <- data.frame(sample_a = unique(res_$sample_a),
                       sample_b = unique(res_$sample_b),
                       barcode_chain_a = unique(res_$barcode_chain_a),
                       barcode_chain_b = unique(res_$barcode_chain_b),
                       sum_a = sum(res_$count_a),
                       sum_b = sum(res_$count_b),
                       total = nrow(res_),
                       significant = sum(res_$significance != "not_significant"),
                       not_significant = sum(res_$significance == "not_significant"),
                       expanded = sum(res_$significance == "B > A"),
                       expanded_freq_a = sum(res_$freq_a[res_$significance == "B > A"]),
                       expanded_freq_b = sum(res_$freq_b[res_$significance == "B > A"]),
                       novel_freq_b = sum(res_$freq_b[res_$significance == "B > A" & res_$count_a == 0]),
                       contracted = sum(res_$significance == "A > B"),
                       contracted_freq_a = sum(res_$freq_a[res_$significance == "A > B"]),
                       contracted_freq_b = sum(res_$freq_b[res_$significance == "A > B"]),
                       lost_freq_a = sum(res_$freq_a[res_$significance == "A > B" & res_$count_b == 0]),
                       novel = sum(res_$significance == "B > A" & res_$count_a == 0),
                       lost = sum(res_$significance == "A > B" & res_$count_b == 0)
    )
    dca_summary <- rbind(dca_summary, res_)
    write.csv(dca_summary, paste0(differential_clone_abundance_results_folder, "/differential_clone_abundance_summary.csv"), row.names = F)
  }
  rm(i)
  return(dca_summary)
}

differential_clone_abundance_summary <- differential_clone_abundance_summarizer()
saveRDS(differential_clone_abundance_summary, paste0(differential_clone_abundance_results_folder, "/differential_clone_abundance_summary.rds"))

print(differential_clone_abundance_summary)

rm(dca_files)




