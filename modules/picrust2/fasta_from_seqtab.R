library(seqinr)
args <- commandArgs(TRUE)

seqtab <- readRDS(args[[1]])

write.fasta(as.list(cols), cols, c(args[[2]]))

