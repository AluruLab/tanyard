#library(lattice)
library(ggplot2)

rle_densityplot <- function(Wd, outfile, sfx){
    #Dirs=dir(Wd, full.names=TRUE)
    y=numeric(0)
    d=Wd
    #for (d in Dirs) {
      fin=paste(paste(d, basename(d), sep="/"), sfx, sep=".")
      x=read.table(fin, head=F)
      y=append(y,x[,3])
      if (max(abs(x[,3])) > 0.5) print(d)
    #}
    # trellis.device(device="pdf", file=outfile)
    # densityplot(y, main=paste("Density of IQR of ", sfx, sep=" "),
    #               xlab = paste(sfx, "values"))
    ysx = summary(y)
    df = data.frame(median = y)
    px = ggplot(df) + geom_density(aes(x = median)) +
         geom_rug(aes(x = median, y = 0), position = position_jitter(height = 0)) +
         xlab(paste(sfx, "; mean = ", as.character(mean(y)), "; sd = ", as.character(sd(y)), sep=""))
    ggsave(filename=outfile, plot=px)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 3) {
    rle_densityplot(args[1], args[2], args[3])
}
