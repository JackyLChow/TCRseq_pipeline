adaptive_results <- read.table("~/Documents/BFX_proj/TCRseq_pipeline/_input/Results_differential-abundance_2023-09-19_03-55/P10-PB1_VS_P10-PB2.differentialAbundance.tsv", header = T, sep = "\t")

# TCRseq_differential_abundance

clone_data_ <- readRDS("~/Documents/BFX_proj/TCRseq_pipeline/_output/clone_data.rds")

a_ <- "P10-PB1.tsv"
b_ <- "P10-PB2.tsv"

counts_A <- clone_data_[clone_data_$file == a_, c("nucleic_acid", "count")]
counts_B <- clone_data_[clone_data_$file == b_, c("nucleic_acid", "count")]
counts <- merge(counts_A, counts_B, by = "nucleic_acid", all = TRUE, suffixes = c("_a","_b"))
counts[is.na(counts)] <- 0

counts6 <- head(counts, 1)

A <- counts$count_a
B <- counts$count_b
sumA <- sum(A)
sumB <- sum(B)

# Morisita index; population similarity, does account for relative abundance
mi <- 2 * sum(A * B/(sumA * sumB)) / (sum((A/sumA)^2) + sum((B / sumB)^2))

# Jaccard similarity; reflects overlap at clone level, does not factor relative abundance
intersect <- length(intersect(counts_A$nucleic_acid, counts_B$nucleic_acid))
union <- length(counts_A$nucleic_acid) + length(counts_B$nucleic_acid) - intersect
js <-  intersect/union
  

# Differential clone abundance
p_value <- c()
for(i in 1:nrow(counts)){
  p_value <- c(p_value,
               binom.test(counts[i, "count_a"], counts[i, "count_a"] + counts[i, "count_b"], p = sumA/(sumA + sumB), alternative = "two.sided")$p.value)
}

differential_abundance <- data.frame(clone = counts[, 1], p_value = p_value)

binom.test(counts6$count_a, counts6$count_a + counts6$count_b, p = sumA/(sumA + sumB), alternative = "two.sided")$p.value

fisher.test(matrix(c(counts6$count_a, sumA, counts6$count_b, sumB), nrow = 2),
            alternative = "two.sided")$p.value
