library(mixtools)

l2variance = function(eqdat){
    gvar = apply(eqdat, 1, var)
    gsm = data.frame("gene" = names(gvar), "var" = gvar)
    gsm$var = log2(gsm$var)  
    gsm
}

select_genes = function(eqdat, nk=3){
    gsm = l2variance(eqdat)
    emrk = normalmixEM(gsm$var, k=nk)
    rjmu = which.min(emrk$mu)
    Midx = which(max.col(emrk$posterior) != rjmu)
    Midx
}

mixt_filter_subset = function(infile, outfile, 
                       start = 500, slab = 500){
    M = read.table(infile, head=T, row.names=1)
    #print()
    nsizes = rev(seq(start, dim(M)[2], by=slab))
    ngenes = rep(0, length(nsizes))
    for(x in (1:length(nsizes))){
        Midx = select_genes(M[, 1:nsizes[x]])
        ngenes[x] = length(Midx)
        M =  M[, 1:nsizes[x]]
        print(c(nsizes[x], ngenes[x]))
    }
    write.table(file=outfile, data.frame(data = nsizes, genes = ngenes),
                sep="\t", quote=F)
}

mixt_filter_sample = function(infile, outfile, 
                              nsample=20,
                              start = 500, slab = 500){
    M = read.table(infile, head=T, row.names=1)
    #print()
    nsizes = rev(seq(start, dim(M)[2], by=slab))
    #ngenes = rep(0, length(nsizes))
    ngenes = matrix(0, nrow = nsizes, ncol = nsample)
    for(x in (1:length(nsizes))){
        for(y in (1:nsample)){
            sx = sample(1:dim(M)[2], nsizes[x])
            Midx = select_genes(M[, sx])
            ngenes[x, y] = length(Midx)
        }
        print(c(nsizes[x], summary(ngenes[x,])))
        M =  M[, 1:nsizes[x]]
    }
    mdf = data.frame(
        cbind(ngenes, t(apply(ngenes, 1 ,summary)))
    )
    mdf$data = nsizes
    write.table(file=outfile, mdf, sep="\t", quote=F)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    mixt_filter_subset(infile=args[1],
                       outfile=args[2])
} else if(length(args) == 3){
    mixt_filter_sample(infile=args[1],
                       outfile=args[2], 
                       nsample=as.numeric(args[3]))
} else {
   print("Usage: Rscript mixt_sampling_filter.R <INFILE> <OUTFILE> [<NSAMPLES>]")
}
