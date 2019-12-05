args <- commandArgs(TRUE)

# LOAD PARAMS
r1 <- c(args[[1]])
r2 <- c(args[[2]])
errF <- readRDS(args[[3]])
errR <- readRDS(args[[4]])
merged <- args[[5]]
minOverlap <-  as.integer(args[[6]])
stats <- c(args[[7]])


library(dada2)

# DEREPLICATE FOR FASTER PROCESSING
derepF <- derepFastq(r1)
derepR <- derepFastq(r2)
derepNum <- length(derepF$map)

# CORE ALGORITHM
ddF <- dada(derepF, err=errF, multithread=TRUE)
ddR <- dada(derepR, err=errR, multithread=TRUE)

# MERGE back
merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap=minOverlap)
mergedNum <- sum(getUniques(merger))

# WRITE MERGED
saveRDS(merger, merged) # CHANGE ME to where you want sequence table saved

A = matrix(c(derepNum, mergedNum), ncol=2)
dimnames(A) = list(c(), c("derep", "merged"))
write.table(A, stats, sep = '\t', row.names = FALSE, quote=FALSE)

