
csv_combine <- function(indir, outfile) {
    files=list.files(path=indir, pattern=".*csv")
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
}
