args <- commandArgs(TRUE)
read_table_loc <- c(args[[1]])
out_loc <- c(args[[2]])
strand <- c(args[[3]])
randomize <- as.logical(args[[4]])
MAX_CONSIST <- as.integer(args[[5]])


library("dada2")
reads <- read.table(file = read_table_loc, sep = '\t', header = TRUE)


if (strand == 'R1'){
    err <- learnErrors(as.character(reads$R1), multithread=TRUE, randomize=randomize, MAX_CONSIST=MAX_CONSIST)
    } else if (strand == 'R2'){
    err <- learnErrors(as.character(reads$R2), multithread=TRUE, randomize=randomize, MAX_CONSIST=MAX_CONSIST)
}

saveRDS(err, out_loc)

