library(stringr)
library(igraph)
source("shinygo_functions.R")

enrich_network = function(genes_lst){
    rgenes = str_c(genes_lst, collapse="\n")
    i = which(names(speciesChoice) == "Arabidopsis thaliana");
    ctd = convertID(rgenes, speciesChoice[[i]])
    # ginfo = geneInfo(ctd, speciesChoice[[i]])
    sgx=significantOverlaps2(rgenes, "Arabidopsis thaliana", "GOBP",
                            0.05, 30)
    if(is.null(sgx)){
        return(NULL)
    }
    g=enrichmentNetwork(sgx,plotting=FALSE)
    g
}

enrich_promoter = function(genes_lst){
    rgenes = str_c(genes_lst, collapse="\n")
    i = which(names(speciesChoice) == "Arabidopsis thaliana");
    ctd = convertID(rgenes, speciesChoice[[i]])
    promoter(ctd, speciesChoice[[i]], "300bp")
}

module.shinygo = function(infile, outfile) {
    rx = read.table(infile, sep="\t", header=T, stringsAsFactors=F)
    mx_mod = max(rx$MODULE_ID)
    rgx = lapply(1:mx_mod, function(mid) rx[rx$MODULE_ID == mid, 'GENE_ID'] )
    print(mx_mod)
    go_net = lapply(rgx, enrich_network)
    go_promoter = lapply(rgx, enrich_promoter)
    ngenes = sapply(rgx, length)
    ncomps = sapply(go_net, function(g) if(is.null(g)){0} else { components(g)$no } )
    ncompsg1 = sapply(go_net, function(g){ if(is.null(g)){0} else { sum(components(g)$csize > 1) }})
    nprom = sapply(go_promoter, function(df){
        if(is.null(df)){
           0
        } else { 
            if(is.null(nrow(df))) 0 else nrow(df)
        } 
    })
    print(length(ngenes))
    print(length(ncomps))
    print(length(ncompsg1))
    print(length(nprom))
    
    #print(ngenes)
    #print(ncomps)
    #print(ncompsg1)
    #print(nprom)
    rsx = data.frame(MODULE_ID=(1:mx_mod), GENES=ngenes, 
                     COMPS=ncomps, COMPGT1=ncompsg1, NPROM=nprom)
            
    write.table(rsx, outfile, sep="\t", )
}

#modfn = "../../data/networks/integ/modules/multilvl/module-assign.tsv"
#print(rsx)

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    module.shinygo(infile=args[1], outfile=args[2])
} else {
   print("Usage: Rscript module_analysis.R <INFILE> <OUTFILE>")
}

dbDisconnect(convert)
