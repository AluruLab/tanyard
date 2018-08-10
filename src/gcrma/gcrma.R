
suppressMessages(library(gcrma))

#outfile = Sys.getenv("R_OUT")
run_gcrma = function(outfile, celfiles){
    D = ReadAffy(filenames=celfiles)
    E = gcrma(D, GSB.adjust=FALSE)
    write.exprs(E, file=outfile)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 1){
    outfile=paste(paste(args[1], basename(args[1]), sep="/"), "gcrma.csv", sep=".")

    run_gcrma(outfile=outfile,
              celfiles=list.celfiles(args[1], full.names=TRUE))
}
