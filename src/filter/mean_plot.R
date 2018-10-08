library(plyr)
library(ggplot2)
library(grid)

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

genes.mean <- function(eqmatrix, log2v=T){
    gsm <- adply(eqmatrix, 1, mean)
    #gsm$xs <- abs(gsm[,5] - gsm[,3])
    rownames(gsm) <- gsm$X1
    colnames(gsm) = c("gene", "mean")
    if(log2v == TRUE){
       gsm = gsm[(gsm$mean > 0), ]
       gsm$mean = log2(gsm$mean)
    }
    gsm
}


mean.deriv.df <- function(gene.summary, binwidth = 0.01,
                        start = 0.01, end = 1.0) {
  mean.breaks <- seq(from = start, to=end, by = binwidth)
  nbreaks <- length(mean.breaks)
  mean.labels <- mean.breaks[2:nbreaks]
  mean.tb <- table(cut(gene.summary$mean,
                      breaks=mean.breaks,
                      labels=mean.labels))
  mean.df <- data.frame(mean = mean.labels, freq = as.vector(mean.tb))
  mean.df
}


genes.mean.plot <- function(gene.summary, start = 0, end = 2,
                           start.limit = 2000, max.limit = 23000){
    ggplot(gene.summary) +
        stat_bin(data = gene.summary, binwidth = 0.01,
                 aes(x = mean, y = 22810 - cumsum(..count..)), geom = "step") +
        xlab("Mean") +
        ylab("# Selected Genes") +
        scale_x_continuous(limits = c(start,end)) +
        scale_y_continuous(limits = c(start.limit, max.limit)) # + theme_bw()
}

mean.bar.plot <- function(mean.df, delta=0, binwidth=0.05,
                         start=0.225, end=1.225, width=0.02){
    tx.df <- if("freq" %in% colnames(mean.df)) { mean.df }
             else { mean.deriv.df(mean.df, binwidth=binwidth, start=start, end=end) }
    tx.df$mean = tx.df$mean + delta
    tx.xlabel = paste("MEAN (± ", binwidth/2, ")", sep ="")
    ggplot(data=tx.df, aes(x=mean, y=freq)) +
        geom_bar(stat='identity', width=width) +
        xlab(tx.xlabel) +
        ylab("No. of Genes")

}

mean.line.plot <- function(mean.df, delta=-0.075, binwidth=0.05,
                         start=0.225, end=1.225){
    #tx.cdf <-  mean.df[mean.df$class==cls.name, ]
    tx.df <- if("freq" %in% colnames(mean.df)) { mean.df }
             else { mean.deriv.df(mean.df, binwidth=binwidth, start=start, end=end) }
    tx.df$mean = tx.df$mean + delta
    tx.xlabel = paste("mean (± ", binwidth/2, ")", sep ="")
    ggplot(data=tx.df, aes(x=mean, y=freq)) +
        geom_point() + geom_line() +
        xlab(tx.xlabel) +
        ylab("No. of Genes")

}

mean.freq.plot <- function(mean.df, binwidth=0.025,
                          start=0.225, end=1.25, width=0.02){
  mean.df = if("freq" %in% colnames(mean.df)) { mean.df }
           else { mean.deriv.df(mean.df, binwidth=binwidth, start=start, end=end) }
  tx.df = mean.df[mean.df$mean <= end, ]
  tx.xlabel = paste("mean (± ", binwidth/2, ")", sep ="")
  ggplot(data=tx.df, aes(x=mean, y=freq)) +
    geom_bar(stat='identity', width=width) +
    xlab(tx.xlabel) +
    ylab("No. of Genes")
}


mean.plot <- function(csvfn, opdf, type = "bar", ... ){
  eqdat = read.table(csvfn, row.names = 1, head = T)
  eqmat = as.matrix(eqdat)
  mean.df = genes.mean(eqmat, FALSE)
  minv = min(mean.df$mean)
  maxv = max(mean.df$mean)
  if(type == "all"){
    p1 = mean.bar.plot(mean.df, start=minv, end=maxv) 
    p2 = mean.line.plot(mean.df, start=minv, end=maxv) 
    p3 = mean.freq.plot(mean.df, start=minv, end=maxv) 
    p4 = genes.mean.plot(mean.df, start=minv, end=maxv) 
    pdf(file=opdf)
    multiplot(p1,p2,p3,p4, cols =2)
    dev.off()
  } else {
    p = if(type == "bar"){ mean.bar.plot(mean.df, start=minv, end=maxv) }
        else if(type == "line") { mean.line.plot(mean.df, start=minv, end=maxv) }
        else if(type == "freq") { mean.freq.plot(mean.df, start=minv, end=maxv) }
        else if(type == "genes") { genes.mean.plot(mean.df, start=minv, end=maxv) }
    ggsave(opdf, plot=p, width = 6, height = 4, units = "in")
  }
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    mean.plot(csvfn=args[1], opdf=args[2])
} else if(length(args) == 3){
    mean.plot(csvfn=args[1], opdf=args[2], type=args[3])
} else {
   print("Usage: Rscript mean_plot.R <INFILE> <OUT_PLOT>")
}
