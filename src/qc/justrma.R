#
# Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
suppressMessages(library(affy))

run_justRMA <- function(celd, outfile){
    celfiles=list.celfiles(path=celd)
    Eset=justRMA(filenames=celfiles, celfile.path=celd)
    D = exprs(Eset)
    D = 2^D
    D = round(D, digits=5)
    write.table(D, file=outfile, sep="\t", quote=F)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2) {
    run_justRMA(args[1], args[2])
}
