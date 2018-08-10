#
# Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
#
run_kstest <- function(infile, outfile){
    D = read.table(infile, head=TRUE, row.names=1)

    m = ncol(D)
    df = m - 2

    M = matrix(nrow=m, ncol=2)

    for (i in 1:m) {
        res = ks.test(D[,i], "pt", df)
        M[i,1] = res$stat
        M[i,2] = 1;
        if (M[i,1] > 0.15) M[i,2] = 0;
    }

    M[,1] = round(M[,1], digits=5)

    write.table(M, outfile, col.names=FALSE, sep="\t", quote=FALSE)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    run_kstest(infile=args[1], outfile=args[2])
}
