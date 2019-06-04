args <- commandArgs(TRUE)
read_table_loc <- c(args[[1]])
out_loc <- c(args[[2]])
strand <- c(args[[3]])

library("dada2")
reads <- read.table(file = read_table_loc, sep = '\t', header = TRUE)


if (strand == 'R1'){
    err <- learnErrors(as.character(reads$R1), multithread=TRUE)
    } else if (strand == 'R2'){
    err <- learnErrors(as.character(reads$R2), multithread=TRUE)
}

saveRDS(err, out_loc)

