library("ath1121501.db")
library("topGO")
library("plyr")

run.top.GO <- function(iqr.vector, go.desc,
                       go.ont = "BP",
                       iqr.limit = 1.0){
    topIQR <- function(x){x > iqr.limit};
    new("topGOdata", description = go.desc, ontology = go.ont,
        allGenes = iqr.vector, geneSel = topIQR,
        nodeSize = 50, annot = annFUN.db, affyLib = "ath1121501.db")
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

top.GO.tests <- function(iqr.vector, go.desc,
                          go.ont = "BP",
                          iqr.limit = 0.8){
    tgo.data <- run.top.GO(iqr.vector, go.desc,
                           go.ont = go.ont, iqr.limit = iqr.limit)
    fisher.results <- run.fisher(tgo.data)
    ks.results <- run.ks(tgo.data)
    elim.ks.results <- run.elim.ks(tgo.data)
    avg.results <- combineResults(fisher.results, ks.results, elim.ks.results)
    list(tgo.data, fisher.results, ks.results, elim.ks.results, avg.results)
}

fisher.top.GO <- function(iqr.vector, go.desc,
                          go.ont = "BP",
                          iqr.limit = 0.8,
                          ntop = 20){
    tgo.data <- run.top.GO(iqr.vector, go.desc,
                           go.ont = go.ont, iqr.limit = iqr.limit)
    fisher.results <- run.fisher(tgo.data)
    ks.results <- run.ks(tgo.data)
    elim.ks.results <- run.elim.ks(tgo.data)
    GenTable(tgo.data, classicFisher = fisher.results,
             classicKS = ks.results, elimKS = elim.ks.results,
             orderBy = "elimKS", ranksOf="classicFisher", topNodes = ntop)
}

common.names <- function(iqr.vector.x, iqr.vector.y, iqr.limit = 0.7){
    union(names(iqr.vector.x[iqr.vector.x > iqr.limit]),
          names(iqr.vector.y[iqr.vector.y > iqr.limit]))
}

common.names.lst <- function(iqr.vec.lst, iqr.limit = 0.7){
    ldply(iqr.vec.lst, function(iqr.vector.x){
        data.frame(probes = as.character(names(
                       iqr.vector.x[iqr.vector.x > iqr.limit])))
    })
}

common.names.lst2 <- function(iqr.vec.lst, iqr.limit = 0.7){
    nlen <- length(iqr.vec.lst)
    ndf <- ldply(1:nlen, function(i){
        iqr.vector.x <- iqr.vec.lst[[i]]
        cnames <- as.character(names(
            iqr.vector.x[iqr.vector.x > iqr.limit]))
        cdf <- data.frame(probes = cnames, iqr = iqr.vector.x[cnames])
        cdf$lst <- i
        cdf
    })
    ndf
    ddply(ndf, .(probes), function(x){
         if(nrow(x) == 1){
             c(iqr = x[1, "iqr"], lst = x[1, "lst"]);
         }
    })
}


fisher.top.common.GO <- function(iqr.vector.x, go.desc.x,
                                 iqr.vector.y, go.desc.y,
                                 go.ont = "BP",
                                 iqr.limit.n = 0.7,
                                 iqr.limit.g = 0.8,
                                 ntop = 20){

    cnames <- common.names(iqr.vector.x, iqr.vector.y,
                           iqr.limit = iqr.limit.n)
    list(fisher.top.GO(iqr.vector.x[cnames], go.desc.x,
                       go.ont = go.ont, iqr.limit = iqr.limit.g,
                       ntop = ntop),
         fisher.top.GO(iqr.vector.y[cnames], go.desc.y,
                       go.ont = go.ont, iqr.limit = iqr.limit.g,
                       ntop = ntop))
}

top.common.lst.GO <- function(iqr.vec.lst, go.desc.x,
                              go.ont = "BP",
                              iqr.limit.n = 0.7,
                              iqr.limit.g = 0.8){
    cnames <- levels(common.names.lst(
        iqr.vec.lst, iqr.limit = iqr.limit.n)$probes)
    run.top.GO(iqr.vec.lst[[1]][cnames], go.desc.x,
               go.ont = go.ont, iqr.limit = iqr.limit.g)
}

top.common.lst.GO.tests <- function(iqr.vec.lst, go.desc.x,
                                     go.ont = "BP",
                                     iqr.limit.n = 1.0,
                                     iqr.limit.g = 1.0){
    cnames <- levels(common.names.lst(
        iqr.vec.lst, iqr.limit = iqr.limit.n)$probes)
    top.GO.tests(iqr.vec.lst[[1]][cnames], go.desc.x,
                 go.ont = go.ont, iqr.limit = iqr.limit.g)
}


fisher.top.common.lst.GO <- function(iqr.vec.lst, go.desc.x,
                                     go.ont = "BP",
                                     iqr.limit.n = 0.7,
                                     iqr.limit.g = 0.8,
                                     ntop = 20){
    cnames <- levels(common.names.lst(
        iqr.vec.lst, iqr.limit = iqr.limit.n)$probes)
    fisher.top.GO(iqr.vec.lst[[1]][cnames], go.desc.x,
                  go.ont = go.ont, iqr.limit = iqr.limit.g,
                  ntop = ntop)
}


top.common.lst2.GO <- function(iqr.vec.lst, go.desc.x,
                               go.ont = "BP",
                               iqr.limit.n = 1.0,
                               iqr.limit.g = 1.0){
    cnames <- as.character(common.names.lst2(
        iqr.vec.lst, iqr.limit = iqr.limit.n)$probes)
    run.top.GO(iqr.vec.lst[[1]][cnames], go.desc.x,
               go.ont = go.ont, iqr.limit = iqr.limit.g)
}


fisher.top.common.lst2.GO <- function(iqr.vec.lst, go.desc.x,
                                      go.ont = "BP",
                                      iqr.limit.n = 1.0,
                                      iqr.limit.g = 1.0,
                                      ntop = 50){
    cnames <- as.character(common.names.lst2(
        iqr.vec.lst, iqr.limit = iqr.limit.n)$probes)
    fisher.top.GO(iqr.vec.lst[[1]][cnames], go.desc.x,
                  go.ont = go.ont, iqr.limit = iqr.limit.g,
                  ntop = ntop)
}

top.common.lst2.GO.tests <- function(iqr.vec.lst, go.desc.x,
                                     go.ont = "BP",
                                     iqr.limit.n = 1.0,
                                     iqr.limit.g = 1.0){
    cnames <- as.character(common.names.lst2(
        iqr.vec.lst, iqr.limit = iqr.limit.n)$probes)
    top.GO.tests(iqr.vec.lst[[1]][cnames], go.desc.x,
                 go.ont = go.ont, iqr.limit = iqr.limit.g)
}
