library(lattice)

rle_densityplot <- function(Wd, outfile){
    Dirs=dir(Wd, full.names=TRUE)
    y=numeric(0)

    for (d in Dirs) {
      fin=paste(paste(d, basename(d), sep="/"), "rle", sep=".")
      x=read.table(fin, head=F)
      y=append(y,x[,3])
      if (max(abs(x[,3])) > 0.5) print(d)
    }
    pdf(outfile, width=20)
    densityplot(y)
    dev.off()
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2) {
    rle_densityplot(args[1], args[2])
}
