args <- commandArgs(TRUE)

fnFs <- c(args[[1]])
fnRs <- c(args[[2]])
filtFs <- c(args[[3]])
filtRs <- c(args[[4]])

truncLenF <- as.integer(args[[5]])
truncLenR <- as.integer(args[[6]])
trimLeftF <- as.integer(args[[7]])
trimLeftR <- as.integer(args[[8]])
maxEEF <- as.numeric(args[[9]])
maxEER <- as.numeric(args[[10]])
truncQ <- as.integer(args[[11]])
maxN <- as.integer(args[[12]])


library(dada2)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncLenF,truncLenR), trimLeft=c(trimLeftF, trimLeftR),
              maxN=maxN, maxEE=c(maxEEF,maxEER), truncQ=truncQ, rm.phix=TRUE, compress=TRUE, minLen = 50)
print(out)