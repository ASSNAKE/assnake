args <- commandArgs(TRUE)

# LOAD PARAMS
read_table_loc <- c(args[[1]])
strand <- c(args[[2]])
derep_out  <- c(args[[3]])

library("dada2")
reads <- read.table(file = read_table_loc, sep = '\t', header = TRUE)

if (strand == 'R1'){
    derep <- derepFastq(as.character(reads$R1))
    saveRDS(derep, derep_out)
} else if (strand == 'R2'){
    derep <- derepFastq(as.character(reads$R2))
    saveRDS(derep, derep_out)
}

