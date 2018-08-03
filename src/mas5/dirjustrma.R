library(affy)

run_dirjustRMA = funciton(Wd){
    setwd(Wd)

    Dirs=dir()

    for (d in Dirs) {
        print(d)
        setwd(d)
        outfile=paste(d, "csv", sep=".")
        Eset=justRMA()
        D=exprs(Eset)
        D=2^D
        D=round(D, digits=5)
        write.table(D, file=outfile, sep="\t", quote=F)
        setwd("..")
    }
}


args = commandArgs(trailingOnly=TRUE)
if(length(args) == 1) {
    cwx = getwd()
    run_dirjustRMA(args[1])
    setwd(cwx)
}
