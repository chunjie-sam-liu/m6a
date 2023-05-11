library(igvR)
igv <- igvR()
setBrowserWindowTitle(igv, "CTCF ChIP-seq")
setGenome(igv, "hg19")
showGenomicRegion(igv, "chr3:128,079,020-128,331,275")



showGenomicRegion(igv, "GATA2")
for(i in 1:4) zoomOut(igv)
