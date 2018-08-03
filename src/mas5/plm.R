library(affyPLM)
library(RColorBrewer)

generate_plots <- function(celfiles, outfile) {
    Cols=brewer.pal(12,"Set3")
    E=ReadAffy(filenames=celfiles)
    Pset=fitPLM(E)
    pdf(outfile, width=20)
    RLE(Pset, main="Relative Log Expression", col=Cols)
    NUSE(Pset, main="Normalized Unscaled Standard Errors", col=Cols)
    dev.off()
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 1){
    generate_plots(celfiles=list.celfiles(path=args[1],
                                          full.names=TRUE),
                   outfile=paste(args[1], "RLE_NUSE.pdf", sep="/"))
}
