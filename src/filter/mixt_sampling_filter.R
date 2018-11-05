library(mixtools)
library(doParallel)
dp.ncores = 1
mixt_filter_sample = function(infile, outfile, 
                              nsample=20,
                              start = 1000,
                              slab = 200,
                              ntop = 10,
                              end = NA) {
 mx.maxit = 4000

 l2variance = function(eqdat){
    gvar = apply(eqdat, 1, var)
    gsm = data.frame("gene" = names(gvar), "var" = gvar)
    gsm$var = log2(gsm$var)  
    gsm
 }

 select_genes = function(eqdat, nk=3){
    gsm = l2variance(eqdat)
    emrk = normalmixEM(gsm$var, k=nk, maxit=4000)
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

 data_dim = function(infile){
    M = read.table(infile, head=T, row.names=1)
    dim(M)
 }

 mixt_sample = function(M, nsz, ntop, nsample, outdir) {
   dM = dim(M)
   topx = matrix(0, nrow = ntop + 1, ncol = dM[1] + 1)
   ngenes = rep(0, nsample)
   for(y in (1:nsample)) {
        sx = sample(1:dM[2], nsz)
        Midx = select_genes(M[, sx])
        ngenes[y] = length(Midx)
        # update best ntop
        topx[ntop + 1, 1:(1 + length(Midx))] = c(length(Midx), Midx)
        rsel = (rank(t(topx[,1]), ties.method="last") > 1)
        topx[1:ntop, ] = topx[rsel, ]
        topx[ntop+1, ] = 0
   }
   print(c(nsz, summary(ngenes)))
   toutf = paste(outdir,
                 paste("mixt-top-", ntop, "-sample-", nsz,  ".csv", sep=""),
                sep="/")
   write.table(file=toutf, topx, sep="\t")
   names(ngenes) = paste("S", 1:length(ngenes))
   c(data = nsz, ngenes, summary(ngenes))
 }


    outdir = dirname(outfile)
    M = read.table(infile, head=T, row.names=1)
    dM = dim(M)
    nsizes = rev(seq(from=start, 
                     to=if(is.na(end)){dM[2]}else{end},
                     by=slab))
    print(nsizes)
    ngenes = matrix(0, nrow = length(nsizes), ncol = nsample)
    # for(x in (1:length(nsizes)))
    cl = makeCluster(dp.ncores)
    registerDoParallel(cl) 
    mdf = foreach(nsz = nsizes, .combine=cbind, .packages='mixtools') %dopar% {
        mixt_sample(M, nsz, ntop, nsample, outdir)
    }
    stopCluster(cl)
    # mdf = data.frame(cbind(ngenes, t(apply(ngenes, 1 ,summary))))
    # mdf$data = nsizes
    write.table(file=outfile, data.frame(t(mdf)), sep="\t", quote=F)
}

args = commandArgs(trailingOnly=TRUE)
ofile = if(length(args) > 2) { paste(args[2], "mixt-samples.csv", sep="/") } else {NA}
if(length(args) == 3){
    mixt_filter_sample(infile=args[1],
                       outfile=ofile, 
                       nsample=as.numeric(args[3]))

} else if(length(args) == 4){
    mixt_filter_sample(infile=args[1],
                       outfile=ofile, 
                       nsample=as.numeric(args[3]),
                       start=as.numeric(args[4]))
} else if(length(args) == 5){
    mixt_filter_sample(infile=args[1],
                       outfile=ofile, 
                       nsample=as.numeric(args[3]),
                       start=as.numeric(args[4]),
                       slab=as.numeric(args[5]))
} else if(length(args) == 6){
    mixt_filter_sample(infile=args[1],
                       outfile=ofile, 
                       nsample=as.numeric(args[3]),
                       start=as.numeric(args[4]),
                       slab=as.numeric(args[5]),
                       ntop=as.numeric(args[6]))
} else if(length(args) == 7){
    mixt_filter_sample(infile=args[1],
                       outfile=ofile, 
                       nsample=as.numeric(args[3]),
                       start=as.numeric(args[4]),
                       slab=as.numeric(args[5]),
                       ntop=as.numeric(args[6]),
                       end=as.numeric(args[7]))
} else {
   print("Usage: Rscript mixt_sampling_filter.R <INFILE> <OUTDIR> [<NSAMPLES> <START> <SLAB> <NTOP> <END>]")
}
