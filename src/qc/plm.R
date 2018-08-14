#
# Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
suppressMessages(library(affyPLM))
suppressMessages(library(RColorBrewer))

generate_plots <- function(celfiles, outfile) {
    Cols=brewer.pal(12,"Set3")
    E=ReadAffy(filenames=celfiles)
    Pset=fitPLM(E)
    pdf(outfile, width=20)
    RLE(Pset, main="Relative Log Expression", col=Cols, las=2, cex.axis=0.4)
    NUSE(Pset, main="Normalized Unscaled Standard Errors", col=Cols, las=2, cex.axis=0.4)
    dev.off()
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 1){
    generate_plots(celfiles=list.celfiles(path=args[1],
                                          full.names=TRUE),
                   outfile=paste(args[1], 
                           paste(basename(args[1]), "RLE_NUSE.pdf", sep="_"), 
                         sep="/"))
}
