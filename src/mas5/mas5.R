
suppressMessages(library(affy))

#outfile = Sys.getenv("R_OUT")
run_mas5 = function(outfile, celfiles){
    D = ReadAffy(filenames=celfiles)
    E = mas5(D, sc=1000)
    write.exprs(E, file=outfile)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 1){
    outfile=paste(paste(args[1], basename(args[1]), sep="/"), "mas5.csv", sep=".")

    run_mas5(outfile=outfile,
             celfiles=list.celfiles(args[1], full.names=TRUE))
}
