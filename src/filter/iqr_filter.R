# infile = Sys.getenv("R_IN")
# outfile = Sys.getenv("R_OUT")
# iqr = as.double(Sys.getenv("R_IQR"))


iqr_filter = function(infile, outfile, iqr){
    print(iqr)
    M=read.table(infile, head=T, row.names=1)
    Iqrs=apply(M, 1, IQR)
    Midx=as.vector(which(Iqrs > iqr, arr.ind=T))

    write.table(file=outfile, M[Midx,], sep="\t", quote=F)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 3){
    iqr_filter(infile=args[1], outfile=args[2], iqr=as.double(args[3]))
} else {
   print("Usage: Rscript iqr_filter.R <INFILE> <OUTFILE> <IQR_VALUE>")
}
