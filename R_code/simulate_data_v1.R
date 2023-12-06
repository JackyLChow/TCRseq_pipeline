# Generate simulated TCRseq data to explore different diversity metrics; 100k total counts, 10k clones

## sample 1: almost perfectly clonal
samp1 <- c(90001, rep(1, 9999))

## sample 2: highly clonal
samp2 <- c(40000, rep(60000/9999, 9999))

## sample 3: oligoclonal
samp3 <- c(40000, 20000, 10000, rep(30000/9997, 9997))

## sample 4: perfect diversity
samp4 <- rep(10, 10000)

