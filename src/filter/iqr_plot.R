suppressMessages(library(plyr))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))

# Multiple plot function
# This function is from cookbook for R
#   http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


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
        xlab("IQR") +
        ylab("# Selected Genes") +
        scale_x_continuous(limits = c(0,2)) +
        scale_y_continuous(limits = c(start.limit, max.limit)) # + theme_bw()
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
    tx.df <- if("freq" %in% colnames(iqr.df)) { iqr.df }
             else { iqr.deriv.df(iqr.df, binwidth=binwidth, start=start, end=end) }
    tx.df$iqr = tx.df$iqr + delta
    tx.xlabel = paste("IQR (± ", binwidth/2, ")", sep ="")
    ggplot(data=tx.df[tx.df$iqr > 0.1, ], aes(x=iqr, y=freq)) +
        geom_bar(stat='identity', width=width) +
        xlab(tx.xlabel) +
        ylab("No. of Genes")

}

iqr.line.plot <- function(iqr.df, delta=-0.075, binwidth=0.05,
                         start=0.225, end=1.225){
    #tx.cdf <-  iqr.df[iqr.df$class==cls.name, ]
    tx.df <- if("freq" %in% colnames(iqr.df)) { iqr.df }
             else { iqr.deriv.df(iqr.df, binwidth=binwidth, start=start, end=end) }
    tx.df$iqr = tx.df$iqr + delta
    tx.xlabel = paste("IQR (± ", binwidth/2, ")", sep ="")
    ggplot(data=tx.df[tx.df$freq > 0, ], aes(x=iqr, y=freq)) +
        geom_point() + geom_line() +
        xlab(tx.xlabel) +
        ylab("No. of Genes")

}

iqr.freq.plot <- function(iqr.df, binwidth=0.025,
                          start=0.225, end=1.25, width=0.02){
  iqr.df = if("freq" %in% colnames(iqr.df)) { iqr.df }
           else { iqr.deriv.df(iqr.df, binwidth=binwidth, start=start, end=end) }
  tx.df = iqr.df[iqr.df$iqr <= end, ]
  tx.xlabel = paste("IQR (± ", binwidth/2, ")", sep ="")
  ggplot(data=tx.df[tx.df$iqr > 0.1, ], aes(x=iqr, y=freq)) +
    geom_bar(stat='identity', width=width) +
    xlab(tx.xlabel) +
    ylab("No. of Genes")
}

iqr.plot <- function(csvfn, opdf, type = "bar", ... ){
  eqdat = read.table(csvfn, row.names = 1, head = T)
  eqmat = as.matrix(eqdat)
  iqr.df = genes.iqr(eqmat)
  #iqr.df = iqr.deriv.df(gsum, ...)
  if(type == "all"){
    p1 = iqr.bar.plot(iqr.df, delta=0) 
    p2 = iqr.line.plot(iqr.df, delta=0) 
    p3 = iqr.freq.plot(iqr.df) 
    p4 = genes.iqr.plot(iqr.df) 
    pdf(file=opdf)
    multiplot(p1,p2,p3,p4, cols =2)
    dev.off()
  } else {
    p = if(type == "bar"){ iqr.bar.plot(iqr.df, delta=0) }
        else if(type == "line") { iqr.line.plot(iqr.df, delta=0) }
        else if(type == "freq") { iqr.freq.plot(iqr.df) }
        else if(type == "genes") { genes.iqr.plot(iqr.df) }
    ggsave(opdf, plot=p, width = 6, height = 4, units = "in")
  }
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    iqr.plot(csvfn=args[1], opdf=args[2])
} else if(length(args) == 3){
    iqr.plot(csvfn=args[1], opdf=args[2], type=args[3])
} else {
   cat("\n")
   cat("  Usage: Rscript iqr_plot.R <INFILE> <OUT_PLOT> [bar/line/freq/genes] \n")
   cat("    <INFILE> : gene expression matrix (GENES X OBS.) with header -- read by read.table \n")
   cat("    <OUT_PLOT>:  output plot file in pdf -- saved with ggplot2::ggsave \n")
   cat("    [bar/line/freq/genes] : Optional : types of plot \n")
   cat("      By default, all plots are generated. \n")
   cat("      bar/line/freq plot the iqr histogram in different ways \n")
   cat("      genes plot is a step plot of the cumulative sum of no. genes \n")
   cat("\n")
}
