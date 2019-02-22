library(mixtools)
library(plyr)

l2variance = function(eqdat){
    eqmat = as.matrix(eqdat) 
    gsm = adply(eqmat, 1, var)
    colnames(gsm) = c("gene", "var")
    gsm$var = log2(gsm$var)  
    gsm
}

mixt_filter = function(infile, outfile, nk=3){
    M = read.table(infile, head=T, row.names=1)
    gsm = l2variance(M)
    #emrk = normalmixEM(gsm$var, k=3)
    emrk = normalmixEM(gsm$var, k=nk, maxit=4000)
    rjmu = which.min(emrk$mu)
    Midx = which(max.col(emrk$posterior) != rjmu)
    print(length(Midx))
    write.table(file=outfile, M[Midx,], sep="\t", quote=F)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    mixt_filter(infile=args[1], outfile=args[2])
} else {
   print("Usage: Rscript mixt_filter.R <INFILE> <OUTFILE> ")
}
