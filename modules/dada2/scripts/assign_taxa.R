args <- commandArgs(TRUE)

# LOAD PARAMS
seqtab <- readRDS(args[[1]])
out <- c(args[[2]])

library("dada2")
taxa <- assignTaxonomy(seqtab, "/data5/bio/databases/dada2/silva_nr_v132_train_set.fa.gz", multithread=20)
saveRDS(taxa, out)

