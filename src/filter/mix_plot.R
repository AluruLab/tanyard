library(mclust)
library(mixtools)
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

sdnorm = function(x, mean=0, sd=1, lambda=1){
    lambda*dnorm(x, mean=mean, sd=sd)
}

l2variance = function(csvfn){
    eqdat = read.table(csvfn, row.names = 1, head = T) 
    eqmat = as.matrix(eqdat) 
    gsm = adply(eqmat, 1, var)
    colnames(gsm) = c("gene", "var")
    gsm$var = log2(gsm$var)  
    gsm
}

mixt.plot = function(gsm, nk){
   emrk = normalmixEM(gsm$var, k=nk)
   p1 = ggplot(data = gsm, aes(x = var)) + 
        geom_histogram(aes(y=..density..), binwidth=0.05,
                       fill="lightgray", color="black")
   fcolors = rep("green", nk)
   fcolors[which.min(emrk$mu)] = "red"
   for (x in (1:nk)) { 
       p1 = p1 + stat_function(fun=sdnorm, args=list(emrk$mu[x], emrk$sigma[x], emrk$lambda[x]), 
                               alpha=0.3, geom="polygon", fill=fcolors[x] ) 
   }
   p1
}


mclust.plot = function(gsm, nk) {
   modx = densityMclust(gsm$var, G=(1:nk)) 
   p1 = ggplot(data = gsm, aes(x = var)) + 
        geom_histogram(aes(y=..density..), binwidth=0.05,
                       fill="lightgray", color="black")
   fcolors = rep("green", nk)
   fcolors[which.min(modx$parameters$mean)] = "red"
   for (x in (1:nk)) { 
       p1 = p1 + stat_function(fun=sdnorm, args=list(modx$parameters$mean[x], 
                                                     sqrt(modx$parameters$variance$sigmasq[x]),
                                                     modx$parameters$pro[x]), 
                               alpha=0.3, geom="polygon", fill=fcolors[x]) 
   }
   p1
} 


args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    gsm = l2variance(csvfn=args[1])
    p1 = mixt.plot(gsm, nk=3)
    p2 = mclust.plot(gsm, nk=3)
    #ggsave(pdfn, plot=p1, width=6, height=6, device="pdf")
    pdf(file=args[2])
    multiplot(p1,p2, cols = 1)
    dev.off()
} else {
   print("Usage: Rscript mix_plot.R <INFILE> <OUT_PLOT>")
}

#mod3 = densityMclust(gsm$var, G=1:3)
#for (x in 1:3) { 
# p2 = p2 + stat_function(fun=sdnorm, args=list(mod3$parameters$mean[x], sqrt(mod3$parameters$variance$sigmasq), mod3$parameters$pro[x]), 
# color="blue") } 
