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

genes.var <- function(eqmatrix){
    gsm <- adply(eqmatrix, 1, var)
    #gsm$xs <- abs(gsm[,5] - gsm[,3])
    rownames(gsm) <- gsm$X1
    colnames(gsm) = c("gene", "var")
    gsm
}


var.deriv.df <- function(gene.summary, binwidth = 0.01,
                        start = 0.01, end = 1.0) {
  var.breaks <- seq(from = start, to=end, by = binwidth)
  nbreaks <- length(var.breaks)
  var.labels <- var.breaks[2:nbreaks]
  var.tb <- table(cut(gene.summary$var,
                      breaks=var.breaks,
                      labels=var.labels))
  var.df <- data.frame(var = var.labels, freq = as.vector(var.tb))
  nrows <- nrow(var.df)
  var.df$csum <- cumsum(var.df$freq)
  var.df$slope <- (c(0, var.df$csum[2:nrows] - var.df$csum[1:(nrows - 1)])) /
      rep(binwidth, nrows)
  var.df$deriv <- - var.df$slope
  var.df
}


genes.var.plot <- function(gene.summary, start.limit = 2000,
                           max.limit = 23000){
    ggplot(gene.summary) +
        stat_bin(data = gene.summary, binwidth = 0.01,
                 aes(x = var, y = 22810 - cumsum(..count..)), geom = "step") +
        xlab("Var") +
        ylab("# Selected Genes") +
        scale_x_continuous(limits = c(0,2)) +
        scale_y_continuous(limits = c(start.limit, max.limit)) # + theme_bw()
}

var.bar.plot <- function(var.df, delta=0, binwidth=0.05,
                         start=0.225, end=1.225, width=0.02){
    tx.df <- if("freq" %in% colnames(var.df)) { var.df }
             else { var.deriv.df(var.df, binwidth=binwidth, start=start, end=end) }
    tx.df$var = tx.df$var + delta
    tx.xlabel = paste("VAR (± ", binwidth/2, ")", sep ="")
    ggplot(data=tx.df[tx.df$var > 0.1, ], aes(x=var, y=freq)) +
        geom_bar(stat='identity', width=width) +
        xlab(tx.xlabel) +
        ylab("No. of Genes")

}

var.line.plot <- function(var.df, delta=-0.075, binwidth=0.05,
                         start=0.225, end=1.225){
    #tx.cdf <-  var.df[var.df$class==cls.name, ]
    tx.df <- if("freq" %in% colnames(var.df)) { var.df }
             else { var.deriv.df(var.df, binwidth=binwidth, start=start, end=end) }
    tx.df$var = tx.df$var + delta
    tx.xlabel = paste("var (± ", binwidth/2, ")", sep ="")
    ggplot(data=tx.df[tx.df$freq > 0, ], aes(x=var, y=freq)) +
        geom_point() + geom_line() +
        xlab(tx.xlabel) +
        ylab("No. of Genes")

}

var.freq.plot <- function(var.df, binwidth=0.025,
                          start=0.225, end=1.25, width=0.02){
  var.df = if("freq" %in% colnames(var.df)) { var.df }
           else { var.deriv.df(var.df, binwidth=binwidth, start=start, end=end) }
  tx.df = var.df[var.df$var <= end, ]
  tx.xlabel = paste("var (± ", binwidth/2, ")", sep ="")
  ggplot(data=tx.df[tx.df$var > 0.27, ], aes(x=var, y=freq)) +
    geom_bar(stat='identity', width=width) +
    xlab(tx.xlabel) +
    ylab("No. of Genes")
}


var.plot <- function(csvfn, opdf, type = "bar", ... ){
  eqdat = read.table(csvfn, row.names = 1, head = T)
  eqmat = as.matrix(eqdat)
  var.df = genes.var(eqmat)
  if(type == "all"){
    p1 = var.bar.plot(var.df) 
    p2 = var.line.plot(var.df) 
    p3 = var.freq.plot(var.df) 
    p4 = genes.var.plot(var.df) 
    pdf(file=opdf)
    multiplot(p1,p2,p3,p4, cols =2)
    dev.off()
  } else {
    p = if(type == "bar"){ var.bar.plot(var.df) }
        else if(type == "line") { var.line.plot(var.df) }
        else if(type == "freq") { var.freq.plot(var.df) }
        else if(type == "genes") { genes.var.plot(var.df) }
    ggsave(opdf, plot=p, width = 6, height = 4, units = "in")
  }
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    var.plot(csvfn=args[1], opdf=args[2])
} else if(length(args) == 3){
    var.plot(csvfn=args[1], opdf=args[2], type=args[3])
} else {
   print("Usage: Rscript var_plot.R <INFILE> <OUT_PLOT>")
}
