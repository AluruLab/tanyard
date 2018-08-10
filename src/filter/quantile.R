library(affy)
library(limma)

# get arguments
#infile = Sys.getenv("R_IN")
#outfile = Sys.getenv("R_OUT")
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
}
