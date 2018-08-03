library(affyPLM)

run_plmstat = function(outfile, celfiles){
    Cel=ReadAffy(filenames=celfiles)
    Pset=fitPLM(Cel)

    # *** RLE ***
    Rle=RLE(Pset, type="stats")

    M=t(Rle)
    M=round(M, digits=5)

    m=nrow(M)
    M=cbind(M, rep(1, m))

    Mean=mean(M[,1])
    Sd=sd(M[,1])

    for (i in 1:m) {
        if (abs(M[i,1]-Mean) > (1.75 * Sd)) M[i,3]=0
        if (abs(M[i,2]) > 0.5) M[i,3]=0
    }

    outRLE=paste(outfile, "rle", sep=".")
    write.table(M, outRLE, col.names=FALSE, sep="\t", quote=FALSE)

    # *** NUSE ***
    Nuse=NUSE(Pset, type="stat")
    M=t(Nuse)
    M=round(M, digits=5)

    m=nrow(M)
    M=cbind(M, rep(1,m))

    Mean=mean(M[,1])
    Sd=sd(M[,1])

    for (i in 1:m) {
        if (abs(M[i,1]-Mean) > (1.75 * Sd)) M[i,3]=0
        if (abs(M[i,2]) > 0.5) M[i,3]=0
    }

    outNUSE=paste(outfile, "nuse", sep=".")
    write.table(M, outNUSE, col.names=FALSE, sep="\t", quote=FALSE)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 1){
    run_plmstat(outfile=args[1],
                celfiles=list.celfiles(full.names=TRUE))
}
if(length(args) > 1){
    run_plmstat(outfile=args[1], celfiles = args[-1])
}
#outfile = Sys.getenv("R_OUT")
