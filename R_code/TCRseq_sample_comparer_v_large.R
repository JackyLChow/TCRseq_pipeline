# adaptive_results <- read.table("~/Documents/BFX_proj/TCRseq_pipeline/_input/Results_differential-abundance_2023-09-19_03-55/P10-PB1_VS_P10-PB2.differentialAbundance.tsv", header = T, sep = "\t")
# adaptive_results$fdr <- p.adjust(adaptive_results$pValue, method = "BH")
# adaptive_results %>%
#   group_by(significance) %>%
#   reframe(min_p = min(pValue), max_p = max(pValue),
#           min_fdr = min(fdr), mad_fdr = max(fdr))

# TCRseq-sample comparer
## Run whole repertoire comparisons

# v_large
## pull individual condensed and downsampled files to run comparisons
differential_clone_abundance_results_folder <- 
  comparison_matrix <- readRDS()

sample_comparer <- function(
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
    
    # sample similarity metrics
    ## Morisita index; population similarity, does account for relative abundance
    mi <- 2 * sum(A * B/(ifelse(sumA == sumB, sumA^2, sumA * sumB))) / (sum((A/sumA)^2) + sum((B / sumB)^2))
    
    ## Jaccard similarity; reflects overlap at clone level, does not factor relative abundance
    intersect <- length(intersect(counts_a[, j_i_clone], counts_b[, j_i_clone]))
    union <- length(counts_a[, j_i_clone]) + length(counts_b[, j_i_clone]) - intersect
    js <-  intersect/union
    
    comparison_results <- rbind(comparison_results,
                                data.frame(sample_a = a_, sample_b = b_,
                                           morisita_index = mi, jaccard_similarity = js,
                                           sum_a = sumA, sum_b = sumB))
    
    cat(paste("completed\n"))
  }
  return(comparison_results)
}

sample_comparison <- sample_comparer(sample_id = "sample_name", j_i_clone = "cdr3_v_j")
print(sample_comparison)
saveRDS(sample_comparison, paste0(output_dir, "comparison_results.rds"))
