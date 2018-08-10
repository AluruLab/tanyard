#infile = Sys.getenv("R_IN")
#outfile = Sys.getenv("R_OUT")

run_submean <- function(infiles, outfiles) {
    Tab=read.table(infile)
    Tab=Tab-apply(Tab, 1, mean)
    Tab=round(Tab, digits=5)
    colnames(Tab)=gsub(".CEL", "", colnames(Tab))
    colnames(Tab)=gsub(".CEL", "", colnames(Tab), ignore.case=FALSE)
    #colnames(Tab)=gsub(".cel", "", colnames(Tab))
    colnames(Tab)=gsub(".gz", "", colnames(Tab))
    write.table(file=outfile, Tab, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    run_submean(infile=args[1], outfile=args[2])
}
