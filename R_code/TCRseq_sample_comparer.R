# adaptive_results <- read.table("~/Documents/BFX_proj/TCRseq_pipeline/_input/Results_differential-abundance_2023-09-19_03-55/P10-PB1_VS_P10-PB2.differentialAbundance.tsv", header = T, sep = "\t")
# adaptive_results$fdr <- p.adjust(adaptive_results$pValue, method = "BH")
# adaptive_results %>%
#   group_by(significance) %>%
#   reframe(min_p = min(pValue), max_p = max(pValue),
#           min_fdr = min(fdr), mad_fdr = max(fdr))

# TCRseq_differential_abundance
clone_data <- readRDS("~/Documents/BFX_proj/TCRseq_pipeline/_output/clone_data.rds")

comparison_matrix <- data.frame(sample_a = "P10-PB1.tsv",
                                sample_b = unique(clone_data$file),
                                id = "nucleic_acid",
                                d_a_method = "binomial")

sample_comparer <- function(
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
    
    counts_a <- clone_dat[clone_dat$file == a_, c(id_, "count")]
    if(any(duplicated(counts_a[, id_]))){
      stop("duplicated ids in counts_A")
    }
    counts_b <- clone_dat[clone_dat$file == b_, c(id_, "count")]
    if(any(duplicated(counts_b[, id_]))){
      stop("duplicated ids in counts_B")
    }
    
    counts <- merge(counts_a, counts_b, by = id_, all = TRUE, suffixes = c("_a","_b"))
    counts[is.na(counts)] <- 0
    
    A <- counts$count_a
    B <- counts$count_b
    sumA <- sum(A)
    sumB <- sum(B)
    
    # sample similarity metrics
    ## Morisita index; population similarity, does account for relative abundance
    mi <- 2 * sum(A * B/(ifelse(sumA == sumB, sumA^2, sumA * sumB))) / (sum((A/sumA)^2) + sum((B / sumB)^2))
    
    ## Jaccard similarity; reflects overlap at clone level, does not factor relative abundance
    intersect <- length(intersect(counts_a$nucleic_acid, counts_b$nucleic_acid))
    union <- length(counts_a$nucleic_acid) + length(counts_b$nucleic_acid) - intersect
    js <-  intersect/union
    
    comparison_results <- rbind(comparison_results,
                                data.frame(sample_a = a_, sample_b = b_, morisita_index = mi, jaccard_similarty = js))
    
    cat(paste("completed\n"))
  }
  return(comparison_results)
}

sample_comparison <- sample_comparer()

