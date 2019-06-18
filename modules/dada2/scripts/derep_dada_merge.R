args <- commandArgs(TRUE)

# LOAD PARAMS
r1 <- c(args[[1]])
r2 <- c(args[[2]])
errF <- readRDS(args[[3]])
errR <- readRDS(args[[4]])
merged <- args[[5]]
minOverlap <-  as.integer(args[[6]])

library(dada2)

# DEREPLICATE FOR FASTER PROCESSING
derepF <- derepFastq(r1)
derepR <- derepFastq(r2)

# CORE ALGORITHM
ddF <- dada(derepF, err=errF, multithread=TRUE)
ddR <- dada(derepR, err=errR, multithread=TRUE)

# MERGE back
merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap=minOverlap)

# WRITE MERGED
saveRDS(merger, merged) # CHANGE ME to where you want sequence table saved

