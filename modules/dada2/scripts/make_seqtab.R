args <- commandArgs(TRUE)

# LOAD PARAMS
mergers <- readRDS(args[[1]])
out <- c(args[[2]])

library("dada2")
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, out)

