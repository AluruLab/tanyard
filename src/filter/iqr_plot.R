
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

genes.iqr.plot <- function(gene.summary, start.limit = 2000,
                           max.limit = 23000){
    ggplot(gene.summary) +
        stat_bin(data = gene.summary, binwidth = 0.01,
                 aes(x = iqr, y = 22810 - cumsum(..count..)), geom = "step") +
        scale_x_continuous(limits = c(0,2)) +
        scale_y_continuous(limits = c(start.limit, max.limit)) +
        theme_bw()
}

make.iqr.plot <- function(iqr.df,
                          font_name = "Calibri",
                          font_size = 18) {
    ggplot(data = iqr.df,aes(x = iqr, y = freq)) +
        scale_x_continuous(limits = c(0.25,1.25)) +
        xlab("IQR") +
        ylab("Genes") +
        geom_point() +
        geom_line() +
        theme(text=element_text(family=font_name, size=font_size)) +
        theme_bw()
}

columns.iqr.bar <- function(iqr.val,
                            font_name = "Calibri",
                            font_size = 18,
                            bwidth = 0.02) {
    ggplot(data = data.frame(iqr = iqr.val), aes(x = iqr)) +
        stat_bin(binwidth = bwidth, geom = "line") +
        xlab("IQR") +
        ylab("Genes") +
        scale_x_continuous(limits = c(0.25,1.25)) +
        theme(text=element_text(family=font_name, size=font_size)) +
        theme_bw()
}

iqr.bar.plot <- function(iqr.df, delta=-0.075, binwidth=0.05,
                         start=0.225, end=1.225, width=0.02){
    #tx.cdf <-  iqr.df[iqr.df$class==cls.name, ]
    tx.df <- iqr.deriv.df(iqr.df, binwidth=binwidth, start=start,
                          end=end)
    tx.df$iqr = tx.df$iqr + delta
    tx.xlabel = paste("IQR (± ", binwidth/2, ")", sep ="")
    ggplot(data=tx.df[tx.df$iqr > 0.38, ], aes(x=iqr, y=freq)) +
        geom_bar(stat='identity', width=width) +
        xlab(tx.xlabel) +
        ylab("No. of Genes")

}

iqr.line.plot <- function(iqr.df, delta=-0.075, binwidth=0.05,
                         start=0.225, end=1.225){
    #tx.cdf <-  iqr.df[iqr.df$class==cls.name, ]
    tx.df <- iqr.deriv.df(iqr.df, binwidth=binwidth, start=start,
                          end=end)
    tx.df$iqr = tx.df$iqr + delta
    tx.xlabel = paste("IQR (± ", binwidth/2, ")", sep ="")
    ggplot(data=tx.df[tx.df$freq > 0, ], aes(x=iqr, y=freq)) +
        geom_point() + geom_line() +
        xlab(tx.xlabel) +
        ylab("No. of Genes")

}

iqr.freq.plot <- function(iqr.df, binwidth=0.025,
                          start=0.225, end=1.25, width=0.02){
  tx.df = iqr.df[iqr.df$iqr <= end, ]
  tx.xlabel = paste("IQR (± ", binwidth/2, ")", sep ="")
  ggplot(data=tx.df[tx.df$iqr > 0.27, ], aes(x=iqr, y=freq)) +
    geom_bar(stat='identity', width=width) +
    xlab(tx.xlabel) +
    ylab("No. of Genes")
}

iqr.plot <- function(csvfn, opdf, type = "bar"){
  eqdat = read.table(csvfn, row.names = 1, head = T)
  eqmat = as.matrix(eqdat)
  gsum = genes.iqr(eqmat)
  if(type == "all"){
    p1 = iqr.bar.plot(gsum) 
    p2 = iqr.line.plot(gsum) 
    p3 = iqr.freq.plot(gsum) 
    p4 = genes.iqr.plot(gsum) 
    pdf(file=opdf)
    multiplot(p1,p2,p3,p4, cols =2)
    dev.off()
  } else {
    p = if(type == "bar"){ iqr.bar.plot(gsum) }
        else if(type == "line") { iqr.line.plot(gsum) }
        else if(type == "freq") { iqr.freq.plot(gsum) }
    ggsave(opdf, plot=p, width = 6, height = 4, units = "in")
  }
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    iqr.plot(csvfn=args[1], opdf=args[2])
} else if(length(args) == 3){
    iqr.plot(csvfn=args[1], opdf=args[2], type=args[3])
} else {
   print("Usage: Rscript iqr_plot.R <INFILE> <OUT_PLOT>")
}
