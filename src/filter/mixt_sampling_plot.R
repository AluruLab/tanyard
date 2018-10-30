library(reshape2)
library(ggplot2)

ndata = c(
"flower"=920,
"leaf"=3575,
"root"=2909,
"rosette"=1311,
"seed"=787,
"seedling1wk"=2834,
"seedling2wk"=1841,
"shoot"=1644,
"wholeplant"=1017
)

mixt_ngenes = c(
"flower"=16165,
"leaf"=16279,
"root"=15689,
"rosette"=15225,
"seed"=15668,
"seedling1wk"=15345,
"seedling2wk"=14114,
"shoot"=15178,
"wholeplant"=15414
)

mixt.table = function(infile, tissue){
    lx = read.table(infile, sep="\t")
    lx
}

mixt.melt.table = function(infile, tissue) {
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
    mlx
}


mixt.plot = function(infile, opdf, tissue) {
    print(sum(ndata))
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
    p1 = ggplot(mlx, aes(x=data, y = ngenes)) + geom_boxplot()
    ggsave(opdf, plot=p1, width = 8, height = 4, units = "in")
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 3){
    mixt.plot(infile=args[1], opdf=args[2], tissue=args[3])
} else {
   print("Usage: Rscript mix_sampling_plot.R <INFILE> <OUT_PLOT> <TISSUE>")
}
