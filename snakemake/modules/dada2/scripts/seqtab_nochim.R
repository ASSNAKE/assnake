args <- commandArgs(TRUE)

# LOAD PARAMS
seqtab <- readRDS(args[[1]])
out <- c(args[[2]])

library("dada2")
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
saveRDS(seqtab_nochim, out)

