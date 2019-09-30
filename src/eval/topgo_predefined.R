library("ath1121501.db")
library("topGO")
library("plyr")
library("stringr")


run.top.GO <- function(iqr.vector, subset.vec, go.desc,
                       go.ont = "BP"){
    boolSelect = names(iqr.vector) %in% subset.vec
    iqrtSelect = iqr.vector >= 0.65
    biqrSelect = boolSelect & iqrtSelect
    print(sum(boolSelect))
    print(sum(iqrtSelect))
    print(quantile(iqr.vector[boolSelect]))
    print(sum(biqrSelect))
    geneList = factor(as.integer(boolSelect))
    geneList = factor(as.integer(biqrSelect))
    names(geneList) = names(iqr.vector)
    new("topGOdata", description = go.desc, ontology = go.ont,
        allGenes = geneList, 
        nodeSize = 50, annot = annFUN.db, affyLib = "ath1121501.db")
}


read.iqr.table = function(iqr.file) {
   tx = read.table(iqr.file, header=T)
   iqrv = tx$IQR; names(iqrv) = tx$PROBE
   iqrv
}



run.fisher <- function(tgo.data){
    runTest(tgo.data, algorithm = "classic", statistic = "fisher")
}

run.ks <- function(tgo.data) {
    runTest(tgo.data, algorithm = "classic", statistic = "ks")
}

run.elim.ks <- function(tgo.data) {
    runTest(tgo.data, algorithm = "elim", statistic = "ks")
}

fisher.top.GO = function(iqr.vec, subset.vec, go.desc,
                         go.ont = "BP", ntop=50){
    tgo.data = run.top.GO(iqr.vec, subset.vec, go.desc)
    fisher.results <- run.fisher(tgo.data)
    ks.results <- run.ks(tgo.data)
    elim.ks.results <- run.elim.ks(tgo.data)
    GenTable(tgo.data, classicFisher = fisher.results,
             classicKS = ks.results, elimKS = elim.ks.results,
             orderBy = "elimKS", ranksOf="classicFisher", topNodes = 20)

}

fisher.diff = function(iqrfile, subfile, outfile, godesc) {
    iqrv = read.iqr.table(iqrfile)
    subgv = read.table(subfile)$V1
    gx.df = fisher.top.GO(iqrv, subgv, godesc)
    write.table(file=outfile, gx.df, sep="\t", quote=F, row.names = F)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 3){
    fisher.diff(iqrfile=args[1], subfile=args[2],
                outfile=args[3], godesc="Predefined")
} else {
 if(length(args) == 4){
    fisher.diff(iqrfile=args[1], subfile=args[2],
                outfile=args[3], godesc=args[4])
 } else {
   print("Usage: Rscript togo_enrichment.R <INFILE1> <SUBFILE> <OUTFILE1> ")
 }
}
