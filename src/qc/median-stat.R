#library(lattice)
library(ggplot2)

rle_densityplot <- function(Wd, outfile, sfx, pth=0){
    Dirs=dir(Wd, full.names=TRUE)
    y=numeric(0)

    for (d in Dirs) {
      fin=paste(paste(d, basename(d), sep="/"), sfx, sep=".")
      if(file.exists(fin)){
        x=read.table(fin, head=F)
        y=append(y,x[,2])
        if (max(abs(x[,2])) > pth) { print(paste("MEDIAN >", pth, basename(d))) }
      } else {
        print(paste("Ignoring...", d))
      }
    }
    #trellis.device(device="pdf", file=outfile)
    #densityplot(y)
    #densityplot(y, main=paste("Density of Meidan of ", sfx, sep=" "),
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
if(length(args) == 4) {
    rle_densityplot(args[1], args[2], args[3], as.numeric(args[4]))
}
