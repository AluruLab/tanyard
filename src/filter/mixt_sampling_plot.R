library(reshape2)
library(ggplot2)
library(grid)
library(dplyr)

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

ndata = c(
"flower"=920,
"leaf"=3575,
"root"=2909,
"rosette"=1311,
"seed"=787,
"seedling1wk"=2834,
"seedling2wk"=1841,
"shoot"=1644,
"wholeplant"=1017,
"complete"=16838
)

mixt_ngenes = c(
"flower"=16335,  # 16165,
"leaf"=16459, #16279,
"root"=15853, #15689,
"rosette"=15387, #15225,
"seed"=15836, #15668,
"seedling1wk"=15503, #15345,
"seedling2wk"=14255, #14114,
"shoot"=15344, #15178,
"wholeplant"=15567, #15414
"complete"=16366
)

mixt.table = function(infile, tissue){
    lx = read.table(infile, sep="\t")
    lx
}

mixt.melt.table = function(infile, tissue,
                           upth = 20000, dwth = 12000) {
    lx = read.table(infile, sep="\t")
    ncols = dim(lx)[2]
    lx = lx[,1:(ncols - 6)] 
    nsamples = dim(lx)[2] - 1
    tx = data.frame(t(c(ndata[tissue],
                    mixt_ngenes[tissue],
                    rep(NA, nsamples - 1))))
    colnames(tx) = c("data", paste("S", 1:nsamples, sep="."))
    lx = rbind(lx,tx)
    lx$data = factor(lx$data)
    mlx = melt(lx, value.name="ngenes")
    mlx = mlx[!is.na(mlx$ngenes), ]
    cat("Removing", sum((mlx$ngenes <= dwth) | (mlx$ngenes >= upth)), "outlier rows\n" )
    mlx = mlx[(mlx$ngenes > dwth) & (mlx$ngenes < upth), ]
    mlx[, c("data", "ngenes")]
}

mixt.melt.summary = function(infile, tissue, 
                             upth = 20000, dwth = 12000) {
    tx = mixt.melt.table(infile, tissue, upth, dwth) %>% 
            group_by(data) %>% 
            summarize(mean = mean(ngenes),
                  median = median(ngenes), 
                  q25 = quantile(ngenes,0.25), 
                  q75=quantile(ngenes,0.75))
    tx$data = as.numeric(levels(tx$data))[ as.numeric(tx$data) ]
    tx
}

mixt.summary.plot = function(infile, tissue, fullpt=TRUE){
    tx = mixt.melt.summary(infile, tissue)
    if(fullpt) {
      ggplot(data=tx[tx$data < max(tx$data),], aes(x=data)) +
        geom_point(aes(y=median, color="median")) +
        geom_smooth(aes(y=median, color="median")) + 
        geom_smooth(aes(y=mean, color="mean")) +
        geom_point(aes(y=mean, color="mean")) + 
        annotate("point", x = ndata[tissue], 
                 y = mixt_ngenes[tissue]) +
        theme(legend.position="bottom")
    } else {
      ggplot(data=tx[tx$data < max(tx$data),], aes(x=data)) +
        geom_point(aes(y=median, color="median")) +
        geom_smooth(aes(y=median, color="median")) + 
        geom_smooth(aes(y=mean, color="mean")) +
        geom_point(aes(y=mean, color="mean")) + 
        theme(legend.position="bottom")
    }
}

mixt.boxplot = function(infile, tissue,
                       upth = 20000, dwth = 12000) {
    mlx = mixt.melt.table(infile, tissue, upth, dwth)
    ggplot(mlx, aes(x=data, y = ngenes)) + 
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 90))
}

mixt.full.plot = function(infile, opdf, tissue, fullpt) {
    p1 = mixt.boxplot(infile, tissue)
    p2 = mixt.summary.plot(infile, tissue, fullpt)
    # tx = mixt.melt.summary(infile, tissue)
    # p2 = ggplot(data=tx[tx$data < max(tx$data),], 
    #             aes(x=data)) +
    #     geom_point(aes(y=median)) +
    #     geom_smooth(aes(y=median)) + 
    #     annotate("point", x = ndata[tissue], 
    #              y = mixt_ngenes[tissue]) +
    #     theme(legend.position="bottom")
    # p3 = ggplot(data=tx[tx$data < max(tx$data),], 
    #             aes(x=data)) +
    #     geom_smooth(aes(y=mean)) +
    #     geom_point(aes(y=mean)) + 
    #     annotate("point", x = ndata[tissue], 
    #              y = mixt_ngenes[tissue]) +
    #     theme(legend.position="bottom")
    pdf(file=opdf, width=6, height=8)
    multiplot(p1, p2, cols=1)
    dev.off()
}

mixt.plot = function(infile, opdf, tissue) {
    print(sum(ndata))
    p1 = mixt.boxplot(infile, tissue)
    ggsave(opdf, plot=p1, width = 8, height = 4, units = "in")
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 3){
    mixt.plot(infile=args[1], opdf=args[2], tissue=args[3])
} else if(length(args) == 4){
    mixt.full.plot(infile=args[1], opdf=args[2], tissue=args[3], fullpt=TRUE)
} else if(length(args) == 5){
    mixt.full.plot(infile=args[1], opdf=args[2], tissue=args[3], fullpt=FALSE)
} else {
   print("Usage: Rscript mix_sampling_plot.R <INFILE> <OUT_PLOT> <TISSUE>")
}
