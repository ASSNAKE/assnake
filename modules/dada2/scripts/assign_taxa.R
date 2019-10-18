args <- commandArgs(TRUE)

# LOAD PARAMS
seqtab <- readRDS(args[[1]])
out <- c(args[[2]])

library("dada2")
taxa <- assignTaxonomy(seqtab, "/data5/bio/databases/dada2/sh_general_release_dynamic_02.02.2019.fa.gz", multithread=20)
saveRDS(taxa, out)

