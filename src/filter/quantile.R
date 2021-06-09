suppressMessages(library(affy))
suppressMessages(library(limma))

# get arguments
quantile_norm <- function(infile, outfile) {
    # get data layout
    D0 = read.table(infile, head=TRUE, row.names=1)
    RNAMES=rownames(D0)
    CNAMES=colnames(D0)
    N = ncol(D0)
    M = nrow(D0)
    write.table(file=outfile, D0, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

    # read and normalize
    D0 = matrix(scan(outfile, n = N * M), M, N, byrow = TRUE)
    D0 = normalizeBetweenArrays(D0, method="quantile")

    D0=round(D0, digits=5)

    # write
    rownames(D0)=RNAMES
    colnames(D0)=CNAMES
    write.table(file=outfile, D0, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
}


args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    quantile_norm(infile=args[1], outfile=args[2])
} else {
   print("Usage: Rscript quantile.R <IN_FILE> <OUT_FILE>")
   print("<IN_FILE>: gene expression matrix (GENES X OBS.) R table fmt. (space seperated w. header)")
   print("<OUT_FILE>: quantile normalized gene expression matrix (GENES X OBS.) R table fmt. (space seperated w. header)")
}
