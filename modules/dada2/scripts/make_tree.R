args <- commandArgs(TRUE)

# LOAD PARAMS
seqtab <- readRDS(args[[1]])
out <- c(args[[2]])
al <- c(args[[3]])

library("dada2")
library("DECIPHER")
library("phangorn")

seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA, processors=12)
saveRDS(alignment, al)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))

saveRDS(fitGTR, out)

