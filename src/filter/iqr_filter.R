# infile = Sys.getenv("R_IN")
# outfile = Sys.getenv("R_OUT")
# iqr = as.double(Sys.getenv("R_IQR"))


function <- iqr_filter(){
    print(iqr)
    M=read.table(infile, head=T, row.names=1)
    Iqrs=apply(M, 1, IQR)
    Midx=as.vector(which(Iqrs > iqr, arr.ind=T))

    write.table(file=outfile, M[Midx,], sep="\t", quote=F)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    iqr_filter(infile=args[1], outfile=args[2], iqrl=as.double(args[3]))
}
