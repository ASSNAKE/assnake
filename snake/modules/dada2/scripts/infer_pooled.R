args <- commandArgs(TRUE)

# LOAD PARAMS
read_table_loc <- c(args[[1]])
err <- readRDS(args[[2]])
strand <- c(args[[3]])
out <- c(args[[4]])
derep_out  <- c(args[[5]])

library("dada2")
reads <- read.table(file = read_table_loc, sep = '\t', header = TRUE)

if (strand == 'R1'){
    derep <- derepFastq(as.character(reads$R1))
    pool <- dada(derep, err=err, multithread=12, pool=TRUE, verbose=0)
    saveRDS(pool, out)
    saveRDS(derep, derep_out)
} else if (strand == 'R2'){
    derep <- derepFastq(as.character(reads$R2))       
    pool <- dada(derep, err=err, multithread=12, pool=TRUE, verbose=0)
    saveRDS(pool, out)
    saveRDS(derep, derep_out)
}

