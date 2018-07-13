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
}

gen_cel_image_df = function(cel_rdf){
    cntr = 0
    lx = purrr::map2(cel_rdf$SampleCEL, cel_rdf$X,
            function(x,y) {
                px1 = str_replace(as.character(x), '.cel.gz', '.pdf');
                px2 = str_replace(px1, '.cel', '.pdf') ;
                cat(as.character(y),",IN,", as.character(x), ",OUT,", px2, "\n");
                gen_cel_image(as.character(x), px2); 
                px2
            })
    data.frame(RowId = cel_rdf$RowId, CEL = cel_rdf$SampleCEL, PDF = unlist(lx))
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