args <- commandArgs(TRUE)

# LOAD PARAMS
rds_fs_loc <- readRDS(args[[1]])
tsv_fs_loc <- c(args[[2]])

library("dada2")

rownames(seqtab) <- gsub('.{3}$', '', sapply(strsplit(rownames(seqtab), ".", fixed = TRUE), `[`, 1))
write.table(seqtab,'./seqtab_nochim.tsv', sep='\t', col.names = NA)


saveRDS(seqtab_nochim, out)
