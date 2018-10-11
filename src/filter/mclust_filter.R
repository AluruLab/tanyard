library(mclust)
library(plyr)

l2variance = function(eqdat){
    eqmat = as.matrix(eqdat) 
    gsm = adply(eqmat, 1, var)
    colnames(gsm) = c("gene", "var")
    gsm$var = log2(gsm$var)  
    gsm
}

mclust_filter = function(infile, outfile){
    M=read.table(infile, head=T, row.names=1)
    gsm = l2variance(M)
    modx = densityMclust(gsm$var, G=(1:3))
    rjmu = which.min(modx$parameters$mean)
    Midx = which(modx$classification != rjmu)
    print(length(Midx))
    write.table(file=outfile, M[Midx,], sep="\t", quote=F)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    mclust_filter(infile=args[1], outfile=args[2])
} else {
   print("Usage: Rscript mclust_filter.R <INFILE> <OUTFILE> ")
}
