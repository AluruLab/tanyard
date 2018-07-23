library("RColorBrewer")
library("affyPLM")
library("purrr")
library("stringr")

gen_cel_image = function(cel_fname, out_fname) {
    rdb_palette <- colorRampPalette(brewer.pal(9, "RdBu"))(32)
    celf <- ReadAffy(filenames = cel_fname)
    pdf(out_fname)
    image(celf[,1], col=rdb_palette)
    dev.off()
    out_fname
}

gen_cel_image_df = function(cel_rdf){
    cntr = 0
    nlen = length(cel_rdf$SampleCEL)
    
    #lx = purrr::map2(cel_rdf$SampleCEL, cel_rdf$X,
    #        function(x,y) {
    lx = rep("NA", nlen)
    lfn = rep("NA", nlen)
    for(idx in (1:nlen)){
            x = cel_rdf$SampleCEL[idx]
            y = cel_rdf$X[idx]
            celfn = as.character(x)
            rowid = as.character(y)
            px1 = str_replace(celfn, '.cel.gz', '.pdf');
            px2 = str_replace(px1, '.cel', '.pdf') ;
            rx = "NA"
            lfn[idx] = px2
            if(file.exists(px2)){
               rx = "PDF FILE EXISTS";
               cat(rowid, ",IN,", celfn, ",OUT,", rx, "\n", sep="");
            } else {
               if(file.exists(celfn)){
                   cat(rowid,",IN,", celfn, ",OUT,"); 
                   rx = tryCatch({gen_cel_image(celfn, px2)}, 
                                  error=function(err) err)
                   if(inherits(rx, "error")) {
                     if(grepl("corrupted", rx)){ 
                       system(paste("rm -f ", celfn));
                       rx =  "CORRUPT CEL FILE";
                     } else {
                       if(grepl("really a CEL file", rx) || grepl("file is truncated", rx)){ 
                         system(paste("rm -f ", celfn));
                         rx = "INVALID CEL FILE";
                       } else{
                         rx = paste("PDF GEN ERROR", as.character(rx));
                       } 
                     } 
                     cat(rx, "\n");
                     lx[idx] = rx
                     next;
                   } else {
                     rx = "CONVERTED"
                   } 
                       
                   cat(rx, "\n");
                   lx[idx] = rx
               } else {
                   rx = "MISSING CEL FILE"
                   cat(rowid,",IN,", celfn, ",OUT,", rx, "\n");
              }
            }
            lx[idx] = rx
    }
    #)
    #data.frame(RowId = cel_rdf$RowId, CEL = cel_rdf$SampleCEL, PDF = unlist(lx))
    data.frame(RowId = cel_rdf$RowId, CEL = cel_rdf$SampleCEL, PDF = lfn, STATUS = lx)
}

gen_cel_image_file = function(in_file, out_file){
    indf = read.csv(in_file, stringsAsFactors=F)
    rdf = gen_cel_image_df(indf)
    write.csv(rdf, out_file)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    cat("Running PDF Generation with", args[1], args[2], "\n")
    gen_cel_image_file(args[1], args[2])
    cat("Completed PDF generation for ", args[1], "\n")
}
