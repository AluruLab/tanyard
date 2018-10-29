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
mclust_ngenes = c(
"flower"=14930,
"leaf"=14897,
"root"=14999,
"rosette"=14559,
"seed"=15380,
"seedling1wk"=15152,
"seedling2wk"=14318,
"shoot"=14800,
"wholeplant"=14371
)

mclust.plot = function(infile, opdf, tissue) {
    print(sum(ndata))
    lx = read.table(infile, sep="\t")
    lx = lx[,c(27, 1:20)] 
    tx = data.frame(t(c(ndata[tissue],
                        mclust_ngenes[tissue],
                        rep(NA, 19))))
    colnames(tx) = c("data", paste("V", 1:20, sep=""))
    lx = rbind(lx,tx)
    lx$data = factor(lx$data)
    mlx = melt(lx, value.name="ngenes")
    mlx = mlx[!is.na(mlx$ngenes), ]
    p1 = ggplot(mlx, aes(x=data, y = ngenes)) + geom_boxplot()
    ggsave(opdf, plot=p1, width = 8, height = 4, units = "in")
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 4){
    mclust.plot(infile=args[1], opdf=args[2], tissue=args[3])
} else {
   print("Usage: Rscript mix_sampling_plot.R <INFILE> <OUT_PLOT> <TISSUE>")
}
