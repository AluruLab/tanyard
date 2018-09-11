#
# Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
suppressMessages(library(affy))

run_dirjustRMA <- function(Wd, pat = NULL){

    Dirs= if (is.null(pat)) { dir(Wd, full.names=TRUE) }
          else { dir(Wd, pattern = pat, full.names=TRUE) }

    for (d in Dirs) {
        outfile=paste(paste(d, basename(d), sep="/"), "csv", sep=".")
        celfiles=list.celfiles(path=d)
        print(paste(d, outfile))
        Eset=justRMA(filenames=celfiles, celfile.path=d)
        D=exprs(Eset)
        D=2^D
        D=round(D, digits=5)
        write.table(D, file=outfile, sep="\t", quote=F)
    }
}


args = commandArgs(trailingOnly=TRUE)
if(length(args) == 1) {
    run_dirjustRMA(args[1])
}
if(length(args) == 2) {
    run_dirjustRMA(args[1], args[2])
}
