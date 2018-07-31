library("RColorBrewer")
library("affyPLM")
library("purrr")
library("stringr")

gen_cel_image_pdf = function(cel_fname, out_fname) {
    rdb_palette <- colorRampPalette(brewer.pal(9, "RdBu"))(32)
    celf <- ReadAffy(filenames = cel_fname)
    pdf(out_fname)
    image(celf[,1], col=rdb_palette)
    dev.off()
    out_fname
}

gen_cel_image_png = function(cel_fname, out_fname) {
    rdb_palette <- colorRampPalette(brewer.pal(9, "RdBu"))(32)
    celf <- ReadAffy(filenames = cel_fname)
    png(out_fname, width=1024, height=1024)
    image(celf[,1], col=rdb_palette)
    dev.off()
    out_fname
}


gen_cel_image = function(cel_fname, out_fname, fmt="pdf") {
    if(fmt == "png"){
       gen_cel_image_png(cel_fname, out_fname)
    } else {
       gen_cel_image_pdf(cel_fname, out_fname)
    } 
}

gen_cel_image_df = function(cel_rdf, fmt="pdf"){
    cntr = 0
    nlen = length(cel_rdf$SampleCEL)
    fmt_ext = ".pdf"
    if(fmt == "png") {
        fmt_ext = ".png"
    }
    
    lx = rep("NA", nlen)
    lfn = rep("NA", nlen)
    for(idx in (1:nlen)){
            x = cel_rdf$SampleCEL[idx]
            y = cel_rdf$X[idx]
            celfn = as.character(x)
            rowid = as.character(y)
            px1 = str_replace(celfn, '.cel.gz', fmt_ext);
            px2 = str_replace(px1, '.cel', fmt_ext) ;
            rx = "NA"
            lfn[idx] = px2
            if(file.exists(px2)){
               rx = "PDF FILE EXISTS";
               cat(rowid, ",IN,", celfn, ",OUT,", rx, "\n", sep="");
            } else {
               if(file.exists(celfn)){
                   cat(rowid,",IN,", celfn, ",OUT,"); 
                   rx = tryCatch({gen_cel_image(celfn, px2, fmt=fmt)}, 
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

    if(fmt == "png"){
       data.frame(RowId = cel_rdf$RowId, CEL = cel_rdf$SampleCEL, PNG = lfn, STATUS = lx)
    } else {
       data.frame(RowId = cel_rdf$RowId, CEL = cel_rdf$SampleCEL, PDF = lfn, STATUS = lx)
    } 
}

gen_cel_image_file = function(in_file, out_file, fmt="pdf"){
    indf = read.csv(in_file, stringsAsFactors=F)
    rdf = gen_cel_image_df(indf, fmt=fmt)
    write.csv(rdf, out_file)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    cat("Running PDF Generation with", args[1], args[2], "\n")
    gen_cel_image_file(args[1], args[2])
    cat("Completed PDF generation for ", args[1], "\n")
}
if(length(args) == 3){
    cat("Running image Generation with", args[1], args[2], args[3], "\n")
    gen_cel_image_file(args[1], args[2], fmt=args[3])
    cat("Completed image generation for ", args[1], " fmt: ", args[3], "\n")
}
