
library(simpleaffy)

run_affyqc <- function(outfile, celfiles) {

    Cel=ReadAffy(filenames = celfiles)
    Qc=qc(Cel)

    M=Qc@average.background
    M=cbind(M, Qc@scale.factors)

    ac3ac5=Qc@qc.probes[,1]/Qc@qc.probes[,2]
    g3g5=Qc@qc.probes[,4]/Qc@qc.probes[,5]

    M=cbind(M, ac3ac5)
    M=cbind(M, g3g5)
    M=round(M, digits=5)
    M=cbind(M, Qc@bioBCalls)

    outqc=paste(outfile, "qc", sep=".")
    write.table(M, outqc, col.names=FALSE, sep="\t", quote=FALSE)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 1){
    run_affyqc(outfile=args[1],
               celfiles=list.celfiles(full.names=TRUE))
}
if(length(args) > 1){
    run_affyqc(outfile=args[1], celfiles = args[-1])
}
#outfile = Sys.getenv("R_OUT")
