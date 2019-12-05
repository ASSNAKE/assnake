args <- commandArgs(TRUE)

# LOAD PARAMS
derepR1 <- readRDS(args[[1]])
derepR2 <- readRDS(args[[2]])
dadaR1 <- readRDS(args[[3]])
dadaR2 <- readRDS(args[[4]])
out <- c(args[[5]])

library("dada2")

mergers <- mergePairs(dadaR1, derepR1, dadaR2, derepR2, verbose=TRUE)

saveRDS(mergers, out)

