#
# Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
suppressMessages(library(affyPLM))
suppressMessages(library(gcrma))

run_gcplmstat = function(outfile, celfiles){
    Cel=ReadAffy(filenames=celfiles)
    CelBG=bg.adjust.gcrma(Cel, GSB.adjust=FALSE)
    Pset=fitPLM(CelBG, background = FALSE)

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

    outRLE=paste(outfile, "gcrle", sep=".")
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

    outNUSE=paste(outfile, "gcnuse", sep=".")
    write.table(M, outNUSE, col.names=FALSE, sep="\t", quote=FALSE)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 1){
    run_gcplmstat(outfile=args[1],
                celfiles=list.celfiles(full.names=TRUE))
}
if(length(args) == 2){
    run_gcplmstat(outfile = args[1],
                celfiles = list.celfiles(args[2], full.names=TRUE))
}
