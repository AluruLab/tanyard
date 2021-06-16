library(plyr)
library(ggplot2)
library(grid)

iqr.deriv.df <- function(gene.summary, binwidth = 0.01,
                        start = 0.01, end = 1.0) {
  iqr.breaks <- seq(from = start, to=end, by = binwidth)
  nbreaks <- length(iqr.breaks)
  iqr.labels <- iqr.breaks[2:nbreaks]
  iqr.tb <- table(cut(gene.summary$iqr,
                      breaks=iqr.breaks,
                      labels=iqr.labels))
  iqr.df <- data.frame(iqr = iqr.labels, freq = as.vector(iqr.tb))
  nrows <- nrow(iqr.df)
  iqr.df$csum <- cumsum(iqr.df$freq)
  iqr.df$slope <- (c(0, iqr.df$csum[2:nrows] - iqr.df$csum[1:(nrows - 1)])) /
      rep(binwidth, nrows)
  iqr.df$deriv <- - iqr.df$slope
  iqr.df
}

genes.iqr <- function(eqmatrix){
    gsm <- adply(eqmatrix, 1, quantile)
    gsm$iqr <- abs(gsm[,5] - gsm[,3])
    rownames(gsm) <- gsm$X1
    gsm
}


iqr.genes <- function(csvfn, outfile) {
  eqdat = read.table(csvfn, row.names = 1, head = T)
  eqmat = as.matrix(eqdat)
  iqr.df = genes.iqr(eqmat)
  iqr.df = iqr.df[, c("X1", "iqr")]
  colnames(iqr.df) = c("PROBE", "IQR")
  write.table(file=outfile, iqr.df, sep="\t", quote=F, row.names = F)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    iqr.genes(csvfn=args[1], outfile=args[2])
} else {
   print("Usage: Rscript iqr_genes.R <IN_FILE> <OUT_FILE>")
   print("<IN_FILE> : gene expression matrix GENES X OBS.")
   print("<OUTFILE> : tab seperated output with PROBE, IQR columns")
}
