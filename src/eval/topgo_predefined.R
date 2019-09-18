library("ath1121501.db")
library("topGO")
library("plyr")
library("stringr")


run.top.GO <- function(iqr.vector, subset.vec, go.desc,
                       go.ont = "BP"){
    geneList = factor(as.integer(subset.vec %in% names(iqr.vector)))
    names(geneList) = names(iqr.vector)
    new("topGOdata", description = go.desc, ontology = go.ont,
        allGenes = iqr.vector, 
        nodeSize = 50, annot = annFUN.db, affyLib = "ath1121501.db")
}


read.iqr.table = function(iqr.file) {
   tx = read.table(iqr.file, header=T)
   iqrv = tx$IQR; names(iqrv) = tx$PROBE
   iqrv
}

fisher.top.GO = function(iqr.vector, subset.vec, go.desc,
                         go.ont = "BP"){
    tgo.data = run.top.GO(iqrv, subgv, go.desc)
    fisher.results <- run.fisher(tgo.data)
    ks.results <- run.ks(tgo.data)
    elim.ks.results <- run.elim.ks(tgo.data)
    GenTable(tgo.data, classicFisher = fisher.results,
             classicKS = ks.results, elimKS = elim.ks.results,
             orderBy = "elimKS", ranksOf="classicFisher", topNodes = ntop)

}

fisher.diff = function(iqrfile, subfile, outfile, godesc) {
    iqrv = read.iqr.table(iqrfile)
    subgv = read.table(iqr.file)$V1
    gx.df = fisher.top.GO(iqrv, subgv, godesc)
    write.table(file=outfile, gx.df, sep="\t", quote=F, row.names = F)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 3){
    fisher.diff(iqrfile=args[1], subfile=args[2],
                outfile=args[3], godesc="Predefined")
} 
else if(length(args) == 4){
    fisher.diff(iqrfile=args[1], subfile=args[2],
                outfile=args[3], godesc=args[4])
} 
else {
   print("Usage: Rscript togo_enrichment.R <INFILE1> <SUBFILE> <OUTFILE1> ")
}
