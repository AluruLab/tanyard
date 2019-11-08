#library(shiny)
library(RSQLite)
library(ggplot2)
#library(grid)
library(gridExtra)

# relative path to data files
datapath = "../../data/shinygo/"   # production server

Min_overlap <- 2
cleanGeneSet = function(x) {
    # remove duplicate; upper case; remove special characters
    x <- unique( toupper( gsub("\n| ","",x) ) )
    x <- x[ which( nchar(x)>1) ]  # genes should have at least two characters
    return(x)
}

# read GMT files, does NO cleaning. Assumes the GMT files are created with cleanGeneSet()
readGMT = function(fileName) {
    x <- scan(fileName, what="", sep="\n")
    x <- strsplit(x, "\t")
    # Extract the first vector element and set it as the list element name
    names(x) <- sapply(x, `[[`, 1)
    x <- lapply(x, `[`, -c(1,2)) # 2nd element is comment, ignored
    x = x[which(sapply(x,length) > 1)]  # gene sets smaller than 1 is ignored!!!
    return(x)
}

sqlite = dbDriver("SQLite")
convert = dbConnect(sqlite,paste0(datapath,"convertIDs.db"),flags=SQLITE_RO)  #read only mode

sqlite  = dbDriver("SQLite")
gmtFiles = list.files(path = paste0(datapath,"pathwayDB"),pattern=".*\\.db")
gmtFiles = paste(datapath,"pathwayDB/",gmtFiles,sep="")
geneInfoFiles = list.files(path = paste0(datapath,"geneInfo"),pattern=".*GeneInfo\\.csv")
geneInfoFiles = paste(datapath,"geneInfo/",geneInfoFiles,sep="")
motifFiles = list.files(path = paste0(datapath,"motif"),pattern=".*\\.db")
motifFiles = paste(datapath,"motif/",motifFiles,sep="")

#STRING10_species = read.csv( paste0(datapath,"data_go/STRING11_species.csv") )

# Create a list for Select Input options
orgInfo <- dbGetQuery(convert, paste("select distinct * from orgInfo " ))
orgInfo <- orgInfo[order(orgInfo$name),]
speciesChoice <- setNames(as.list( orgInfo$id ), orgInfo$name2 )
# add a defult element to list    # new element name       value
speciesChoice <- append( setNames( "BestMatch","Best matching species"), speciesChoice  )
# move one element to the 2nd place
move2 <- function(i) {
    c(speciesChoice[1],speciesChoice[i],speciesChoice[-c(1,i)])
}
i = which(names(speciesChoice) == "Glycine max"); speciesChoice <- move2(i)
i = which(names(speciesChoice) == "Zea mays"); speciesChoice <- move2(i)
i = which(names(speciesChoice) == "Arabidopsis thaliana"); speciesChoice <- move2(i)
i = which(names(speciesChoice) == "Saccharomyces cerevisiae"); speciesChoice <- move2(i)
i = which(names(speciesChoice) == "Caenorhabditis elegans"); speciesChoice <- move2(i)
i = which(names(speciesChoice) == "Zebrafish" ); speciesChoice <- move2(i)
i = which(names(speciesChoice) == "Cow" ); speciesChoice <- move2(i)
i = which(names(speciesChoice) == "Rat" ); speciesChoice <- move2(i)
i = which(names(speciesChoice) == "Mouse"); speciesChoice <- move2(i)
i = which(names(speciesChoice) == "Human"); speciesChoice <- move2(i)

GO_levels = dbGetQuery(convert, "select distinct id,level from GO
                                 WHERE GO = 'biological_process'" )
level2Terms = GO_levels[which(GO_levels$level %in% c(2,3))  ,1]  # level 2 and 3

idIndex = dbGetQuery(convert, paste("select distinct * from idIndex "))

quotes = dbGetQuery(convert, " select * from quotes")
quotes = paste0("\"",quotes$quotes,"\"", " -- ",quotes$author,".       ")

proper = function(x) {
    paste0(toupper(substr(x, 1, 1)), substring(x, 2))
}

extract1 = function (x) {
  words <- unlist (strsplit(x,"_"))
  if(length(words)  <=4) return(gsub("_"," ",x)) else {
  words <- words[-c(1:4)]
  return(proper(paste(words,collapse = " ")))}
}

#find idType based on index
findIDtypeById = function(x){ # find
    return(idIndex$idType[ as.numeric(x)])
}

findSpeciesById = function (speciesID){ # find species name use id
    return(orgInfo[which(orgInfo$id == speciesID),] )
}

# just return name
findSpeciesByIdName = function (speciesID){ # find species name use id
    return(orgInfo[which(orgInfo$id == speciesID),3] )
}

#Homo sapies --> hsapiens
shortSpeciesNames = function(tem){
    tem2 = strsplit(as.character(tem)," ") 	   
    return(tolower(paste0(substr(tem2[[1]][1],1,1), tem2[[1]][2] )))
}

# convert sorted species:idType combs into a list for repopulate species choice
matchedSpeciesInfo = function (x) {
    a<- c()
    for( i in 1:length(x)) {
        a = c(a,paste( gsub("genes.*","",findSpeciesByIdName( as.numeric(gsub(" .*","",names(x[i])) ))), " (",
                    x[i]," mapped from ",findIDtypeById( gsub(".* ","",names(x[i]) ) ),")",sep="")
        )
    }
    return(a)
}

# convert gene IDs to ensembl gene ids and find species
convertID <- function (query, selectOrg) {
    query <- gsub("\"|\'","",query) 
    # remove " in gene ids, mess up SQL query				
    # remove ' in gene ids				
    # |\\.[0-9] remove anything after A35244.1 -> A35244  
    #  some gene ids are like Glyma.01G002100

    querySet = cleanGeneSet(unlist(strsplit(toupper(query),'\t| |\n|\\,')))
 	
    result = dbGetQuery( convert,
        paste(" select distinct id,ens,species from mapping where id IN ('",
              paste(querySet,collapse="', '"),   "')",sep=""))
    if(dim(result)[1] == 0){
        return(NULL)
    }
  
    if(selectOrg == speciesChoice[[1]]) {
        comb = paste(result$species, result$idType)
        sortedCounts = sort(table(comb), decreasing=T)
        recognized = names(sortedCounts[1]  )
        result <- result[which(comb == recognized )  , ]

        speciesMatched = sortedCounts
        names(speciesMatched) = sapply(as.numeric(gsub(" .*","",names(sortedCounts))),
                                       findSpeciesByIdName)
        speciesMatched <- as.data.frame( speciesMatched )
        if(length(sortedCounts) == 1) { # if only  one species matched
            speciesMatched[1,1] <-paste( rownames(speciesMatched), "(",speciesMatched[1,1],")",sep="")
        } else {# if more than one species matched
            speciesMatched[,1] <- as.character(speciesMatched[,1])
            speciesMatched[,1] <- paste(speciesMatched[,1]," (",speciesMatched[,2], ")", sep="")
            speciesMatched[1,1] <- paste(speciesMatched[1,1],
                "   ***Used in mapping***  To change, select from above and resubmit query.")
            speciesMatched <- as.data.frame(speciesMatched[,1])
        }

    } else { # if species is selected
        result <- result[which(result$species == selectOrg ) ,]
        if( dim(result)[1] == 0  ) return(NULL) #stop("ID not recognized!")
        speciesMatched <- as.data.frame(paste("Using selected species ",
                                              findSpeciesByIdName(selectOrg)))
    }
    result <- result[which(!duplicated(result[,1])), ] # remove duplicates in query gene ids 
    result <- result[which(!duplicated(result[,2])), ] # remove duplicates in ensembl_gene_id  
    colnames(speciesMatched) = c("Matched Species (genes)")
    conversionTable <- result[,1:2]
    colnames(conversionTable) = c("User_input","ensembl_gene_id")
    conversionTable$Species = sapply(result[,3], findSpeciesByIdName)
    return(list(originalIDs = querySet,
                IDs=unique( result[,2]),
                species = findSpeciesById(result$species[1]),
                #idType = findIDtypeById(result$idType[1] ),
                speciesMatched = speciesMatched,
                conversionTable = conversionTable))
}

geneInfo <- function (converted, selectOrg){
    if(is.null(converted)){ # no ID
        return(as.data.frame("ID not recognized!"))
    }
    querySet <- converted$IDs
    if(length(querySet) == 0){
        return(as.data.frame("ID not recognized!") )
    }
    ix = grep(converted$species[1,1],geneInfoFiles)
    if (length(ix) == 0 ) {
        return(as.data.frame("No matching gene info file found") )
    } else {
        # If selected species is not the default "bestMatch",
        # use that species directly
        if(selectOrg != speciesChoice[[1]]) {
            ix = grep(findSpeciesById(selectOrg)[1,1], geneInfoFiles)
        }
        if(length(ix) == 1){  # if only one file           #WBGene0000001 some ensembl gene ids in lower case
            x = read.csv(as.character(geneInfoFiles[ix]) ); x[,1]= toupper(x[,1])
        }
        else{ # read in the chosen file
            return(as.data.frame("Multiple geneInfo file found!"))
        }

        Set = match(x$ensembl_gene_id, querySet)
        Set[which(is.na(Set))]="Genome"
        Set[which(Set!="Genome")] ="List"
        # x = cbind(x,Set) } # just for debuging
        return( cbind(x,Set) )
    }
}

promoter <- function (converted, selectOrg, radio){
  idNotRecognized = as.data.frame("ID not recognized!")
  
  if(is.null(converted) ) 
    return(idNotRecognized) # no ID
  
  querySet <- converted$IDs;
  
  if(length(querySet) == 0) 
    return(idNotRecognized )
  ix = grep(converted$species[1,1],motifFiles)

  # If selected species is not the default "bestMatch", use that species directly
  if(selectOrg != speciesChoice[[1]]) {
    ix = grep(findSpeciesById(selectOrg)[1,1], motifFiles )
  }
    
  ix1 =grep(as.character(radio),motifFiles[ix]) # match 300bp or 600bp
  if(length(ix1) >0) ix = ix[ix1]   # if 600 is not found, use 300bp
  if (length(ix) == 0 ) {
    return(as.data.frame("No matching motif file found") )
  } 
  else {
    
    if(length(ix) > 1)  # if only one file
      return(as.data.frame("Multiple geneInfo file found!") )
  
    motifs <- dbConnect(sqlite,motifFiles[ix]) # makes a new file
  
    sqlQuery = paste( " select * from scores where row_names IN ('", paste(querySet,collapse="', '"),"')" ,sep="")
    result <- dbGetQuery( motifs, sqlQuery  )
    if( dim(result)[1] ==0) {
        dbDisconnect(motifs)
        return(list( x=as.data.frame("No matching species or gene ID file!" )) )
    }
    row.names(result) <- result$row_names; result <- result[,-1]
    TFstat <- as.data.frame( cbind(apply(result,2,mean),apply(result,2,sd) ) )
    colnames(TFstat) = c("scoreMean1","scoreSD1" )
    rownames(TFstat) = toupper( colnames(result) )
  
    TFs <- dbGetQuery(motifs, "select ID,TF_Name,Family_Name,DBID,Motif_ID,coreMotif,memo,nGenes,scoreSD,scoreMean from  TF_Information ")
    dbDisconnect(motifs)
    TFs$ID <- toupper(TFs$ID)
  
    TFs <- merge(TFs, TFstat, by.x = 'ID', by.y='row.names')
    TFs <- TFs[!is.na(TFs$scoreSD) ,]  #some TFs return NA -Inf
    n1 = dim(result)[1] # number of genes in query set
    TFs$scoreMean2 <- (TFs$scoreMean * TFs$nGenes - TFs$scoreMean1 *n1)/(TFs$nGenes - n1)
    #SD2 needs to be adjusted too, but ignored for now. use overall SD2
    # t test unequal variance statistic
    TFs$t <- (TFs$scoreMean1-TFs$scoreMean2)/ sqrt( TFs$scoreSD1^2/n1 + TFs$scoreSD^2/TFs$nGenes   )
    # degree of freedom
    TFs$df <- ( TFs$scoreSD1^2/n1 + TFs$scoreSD^2/TFs$nGenes)^2 / ((TFs$scoreSD1^2/n1)^2/(n1-1) +   (TFs$scoreSD^2/TFs$nGenes)^2/(TFs$nGenes-1))
    TFs$pVal =1-pt(TFs$t,df = TFs$df)  # t distribution
    TFs$FDR = p.adjust(TFs$pVal,method="fdr")
    TFs <- TFs[order(TFs$pVal) ,]
    TFs$scoreDiff = round(TFs$scoreMean1 - TFs$scoreMean2,0)
    #TFs <- TFs[order(-TFs$scoreDiff) ,]
  
    # does this transcription factor gene in this cluster?
    ix <- match(toupper( TFs$DBID), querySet) # assuming the DBID column in cisbp are ensembl gene ids
    TFs$note = ""
    if(sum(!is.na(ix)) >0) {
      TFs$note[which(!is.na(ix))] <- "* Query Gene"
    }
    TFs <- subset(TFs, FDR<0.25, select=c(coreMotif,TF_Name,Family_Name, pVal,FDR,scoreDiff, note ) )
    colnames(TFs) =c("Enriched motif in promoter", "TF","TF family","P val.","FDR","Score","Note"   )
    if(dim(TFs)[1] >30 ) 
      TFs <- TFs[1:30,]
    if(dim(TFs)[1] ==0) 
      return(as.data.frame("No significant TF binding motif detected.") ) else
    return( TFs )
   }
}

# Main function. Find a query set of genes enriched with functional category
FindOverlap <- function (converted, gInfo, GO, selectOrg, minFDR, input_maxTerms) {
    idNotRecognized = list(x=as.data.frame("ID not recognized!"),
                            groupings= as.data.frame("ID not recognized!")  )
    if(is.null(converted) ) { # no ID
        return(idNotRecognized)
    }
    querySet <- converted$IDs;
    if(length(querySet) == 0) return(idNotRecognized )

    ix = grep(converted$species[1,1],gmtFiles)
    totalGenes <- converted$species[1,7]

    errorMessage = list(x=as.data.frame("Annotation file cannot be found"),
                        groupings= as.data.frame("Annotation file cannot be found")  )
    if (length(ix) == 0) {
        return( errorMessage )
    }

    # If selected species is not the default "bestMatch", use that species directly
    if(selectOrg != speciesChoice[[1]]) {
        ix = grep(findSpeciesById(selectOrg)[1,1], gmtFiles)
        if (length(ix) == 0 ) {
            return(idNotRecognized)
        }
        totalGenes <- orgInfo[which(orgInfo$id == as.numeric(selectOrg)),7]
    }
    pathway <- dbConnect(sqlite,gmtFiles[ix])

    # Generate a list of geneset categories such as "GOBP", "KEGG" from file
    geneSetCategory <-  dbGetQuery(pathway, "select distinct * from categories " )
    geneSetCategory  <- geneSetCategory[,1]
    categoryChoices <- setNames(as.list( geneSetCategory ), geneSetCategory )
    categoryChoices <- append( setNames( "All","All available gene sets"), categoryChoices  )
    #change GOBO to the full description for display
    names(categoryChoices)[ match("GOBP",categoryChoices)  ] <- "GO Biological Process"
    names(categoryChoices)[ match("GOCC",categoryChoices)  ] <- "GO Cellular Component"
    names(categoryChoices)[ match("GOMF",categoryChoices)  ] <- "GO Molecular Function"

    sqlQuery = paste(" select distinct gene,pathwayID from pathway where gene IN ('",
                     paste(querySet,collapse="', '"),"')" ,sep="")

    #cat(paste0("HH",GO,"HH") )

    if( GO != "All")
        sqlQuery = paste0(sqlQuery, " AND category ='",GO,"'")
    result <- dbGetQuery( pathway, sqlQuery  )
    if( dim(result)[1] ==0) {
        dbDisconnect(pathway)
        return(list( x=as.data.frame("No matching species or gene ID file!" )) )
    }

    # given a pathway id, it finds the overlapped genes, symbol preferred
    sharedGenesPrefered <- function(pathwayID) {
        tem <- result[which(result[,2]== pathwayID ),1]
        ix = match(tem, converted$conversionTable$ensembl_gene_id) # convert back to original
        tem2 <- unique( converted$conversionTable$User_input[ix] )
        if(length(unique(gInfo$symbol) )/dim(gInfo)[1] >.7  ) # if 70% genes has symbol in geneInfo
	    { 
            ix = match(tem, gInfo$ensembl_gene_id);
	        tem2 <- unique( gInfo$symbol[ix] )
        }
        return( paste( tem2 ,collapse=" ",sep="") )
    }

    x0 = table(result$pathwayID)
    x0 = as.data.frame( x0[which(x0>=Min_overlap)] )# remove low overlaps
    errorMessage = list(x=as.data.frame("Too few genes."),
                        groupings= as.data.frame("Too few genes.")  )
    if(dim(x0)[1] <= 5 ) return(errorMessage) # no data
    colnames(x0)=c("pathwayID","overlap")
    pathwayInfo <- dbGetQuery(pathway,
        paste( " select distinct id,n,Description from pathwayInfo where id IN ('",
		        paste(x0$pathwayID,collapse="', '"),   "') ",sep="") )
    x = merge(x0,pathwayInfo, by.x='pathwayID', by.y='id')

    x$Pval=phyper(x$overlap-1,length(querySet),totalGenes - length(querySet),as.numeric(x$n), lower.tail=FALSE );
    x$FDR = p.adjust(x$Pval,method="fdr")
    x <- x[ order( x$FDR)  ,]  # sort according to FDR

    # Gene groups for high level GOBP terms
    groups <- dbGetQuery( pathway,paste( " select distinct id, description from pathwayInfo
                         where golevel IN ( '2','3') ",sep="") )

    ix = match(groups$id, x0$pathwayID)
    if(length(groups)>0 && length(ix) >0) {
        groupings = as.data.frame("No grouping.")
    }
    groups$ngenes <- x0$overlap[ix]
    groups <- groups[which(!is.na(ix) ),]
    groups <- groups[order(-groups$ngenes),]
    if(max(groups$ngenes)<=2) {
        groups = as.data.frame("Too few genes")
    } else {
		groupings = subset(groups,ngenes>2) # at least 10 genes
		if(dim(groups)[1] > 100){
            groups <- groups[1:100,]
        }
        groups = cbind(groups, sapply( groups$id, sharedGenesPrefered ) )
        groups = groups[,-1]
        groups = groups[,c(2,1,3)]
        colnames(groups) = c("N","High level GO category", "Genes")
	}

    if(min(x$FDR) > minFDR){
        x=as.data.frame("No significant enrichment found!")
    } else {
        x <- x[which(x$FDR < minFDR),]
        if(dim(x)[1] > as.integer(input_maxTerms) ){
            x = x[ 1:as.integer(input_maxTerms), ]
        }
        x = cbind(x,sapply( x$pathwayID, sharedGenesPrefered ) )
        colnames(x)[7]= "Genes"
        x <- subset(x,select = c(FDR,overlap,n,description,Genes) )
        colnames(x) = c("Enrichment FDR", "Genes in list", "Total genes","Functional Category","Genes"  )
    }

    dbDisconnect(pathway)
    return(list( x=x, groupings = groups, categoryChoices = categoryChoices ) )
}

significantOverlaps = function(ingenes, orgName,  GOID = "GOBP", minFDR = 0.05, maxTerms = 30){
    i = which(names(speciesChoice) == orgName);
    ctd = convertID(ingenes, speciesChoice[[i]])
    ginfo = geneInfo(ctd, speciesChoice[[i]])
    FindOverlap(ctd, ginfo, GOID, speciesChoice[[i]], minFDR, maxTerms)
}

significantOverlaps2 = function(ingenes, orgName,  GOID = "GOBP", minFDR = 0.05, maxTerms = 30){
    tem = significantOverlaps(ingenes, orgName, GOID, minFDR, maxTerms)
    if(dim(tem$x)[2] ==1 ) return(NULL)
	tem <- tem$x;
    colnames(tem)= c("adj.Pval","nGenesList","nGenesCategor","Pathways","Genes")
    tem$Direction ="Diff"	
    tem
}

# numChar=100 maximum number of characters
# n=200  maximum number of nodes
# degree.cutoff = 0    Remove node if less connected
#from PPInfer
enrich.net2 <-  function (x, gene.set, node.id, node.name = node.id, pvalue, 
    n = 50, numChar = NULL, pvalue.cutoff = 0.05, edge.cutoff = 0.05, 
    degree.cutoff = 0, edge.width = function(x) {
        5 * x^2
    }, node.size = function(x) {
        2.5 * log10(x)
    }, group = FALSE, group.color = c("green","red" ), group.shape = c("circle", 
        "square"), legend.parameter = list("topright"), show.legend = TRUE,
    plotting=TRUE, layoutButton = 0,
    ...) 
{
	library(igraph)
	set.seed(layoutButton)
    x <- data.frame(x, group)
    colnames(x)[length(colnames(x))] <- "Group"
    x <- x[as.numeric( x[, pvalue]) < pvalue.cutoff, ]
    x <- x[order(x[, pvalue]), ]
    n <- min(nrow(x), n)
    if (n == 0) {
        stop("no enriched term found...")
    }
    x <- x[1:n, ]
    index <- match(x[, node.id], names(gene.set))
    geneSets <- list()
    for (i in 1:n) {
        geneSets[[i]] <- gene.set[[index[i]]]
    }
    names(geneSets) <- x[, node.name]
    if (is.null(numChar)) {
        numChar <- max(nchar(as.character(x[, node.name])))
    }
    else {
        if (length(unique(substr(x[, node.name], 1, numChar))) < 
            nrow(x)) {
            numChar <- max(nchar(as.character(x[, node.name])))
            message("Note : numChar is too small.", "\n")
        }
    }
    x[, node.name] <- paste(substr(x[, node.name], 1, numChar), 
        ifelse(nchar(as.character(x[, node.name])) > numChar, 
            "...", ""), sep = "")
    w <- matrix(NA, nrow = n, ncol = n)

    for (i in 1:n) {
        for (j in i:n) {
            u <- unlist(geneSets[i])
            v <- unlist(geneSets[j])
            w[i, j] = length(intersect(u, v))/length(unique(c(u, 
                v)))
        }
    }
    list.edges <- stack(data.frame(w))
    list.edges <- cbind(list.edges[, 1], rep(x[, node.name], 
        n), rep(x[, node.name], each = n))
    list.edges <- list.edges[list.edges[, 2] != list.edges[,3], ]
    list.edges <- list.edges[!is.na(list.edges[, 1]), ]
    if(is.null(list.edges)) return (NULL)
    if(is.null(dim(list.edges))) return (NULL)
    #print(dim(list.edges))
    g <- graph.data.frame(list.edges[, -1], directed = FALSE)
    E(g)$width = edge.width(as.numeric(list.edges[, 1]))
    V(g)$size <- node.size(lengths(geneSets))
    g <- delete.edges(g, E(g)[as.numeric(list.edges[, 1]) < edge.cutoff])
    index.deg <- igraph::degree(g) >= degree.cutoff
    g <- delete.vertices(g, V(g)[!index.deg])
    x <- x[index.deg, ]
    index <- index[index.deg]
    if (length(V(g)) == 0) {
        # stop("no categories greater than degree.cutoff...")
        return (NULL)
    }
    n <- min(nrow(x), n)
    x <- x[1:n, ]
    group.level <- sort(unique(group))
    pvalues <- x[, pvalue]
    for (i in 1:length(group.level)) {
        index <- x[, "Group"] == group.level[i]
        V(g)$shape[index] <- group.shape[i]
        group.pvalues <- pvalues[index]
        if (length(group.pvalues) > 0) {
            if (max(group.pvalues) == min(group.pvalues)) {
                V(g)$color[index] <- adjustcolor(group.color[i], 
                  alpha.f = 0.5)
            }
            else {
                V(g)$color[index] <- sapply(1 - (group.pvalues - 
                  min(group.pvalues))/(max(group.pvalues) - min(group.pvalues)), 
                  function(x) {
                    adjustcolor(group.color[i], alpha.f = x)
                  })
            }
        }
    }
	if(plotting) { 
		plot(g,, vertex.label.dist=0.8, ...)
		if (show.legend) {
			legend.parameter$legend <- group.level
			legend.parameter$text.col <- group.color
			legend.parameter$bty <- "n"	
			do.call(legend, legend.parameter)
		}
    }
    return(g)
}


enrichmentNetwork <- function(enrichedTerms,layoutButton=0,plotting=TRUE){
	geneLists = lapply(enrichedTerms$Genes, function(x) unlist( strsplit(as.character(x)," " )   ) )
	names(geneLists)= enrichedTerms$Pathways
	enrichedTerms$Direction = gsub(" .*","",enrichedTerms$Direction )

	g <- enrich.net2(enrichedTerms, geneLists, node.id = "Pathways", numChar = 100, 
	   pvalue = "adj.Pval", edge.cutoff = 0.2, pvalue.cutoff = 1, degree.cutoff = 0,
	   n = 200, group = enrichedTerms$Direction, vertex.label.cex = 1, vertex.label.color = "black",
       show.legend=FALSE, layoutButton=layoutButton, plotting=plotting)
    g
}
