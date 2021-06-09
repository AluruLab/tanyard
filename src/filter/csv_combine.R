
csv_combine <- function(indir, outfile) {
    files=list.files(path=indir, pattern=".*csv", full.names=T)
    tab=read.table(files[1])
    for (i in files[-1]) {
        print(i)
        x=read.table(i)
        tab=cbind(tab,x)
    }
    write.table(tab, file=outfile, sep="\t", quote=FALSE)
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
    csv_combine(indir=args[1], outfile=args[2])
} else {
    print("Usage: Rscript csv_combine.R <IN_FILE> <OUT_FILE>")
    print("<IN_FILE>: csv files generated from processing an experiment")
    print("<OUT_FILE>: combined gene expression matrix (GENES X OBS.) in R table fmt.")
}
