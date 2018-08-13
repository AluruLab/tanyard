
#
# Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>

compute_summary <- function(d){
    fpx = paste(d, basename(d), sep="/")
    Rle=read.table(paste(fpx, "rle", sep="."), head=F)
    Nuse=read.table(paste(fpx, "nuse", sep="."), head=F)
    Qc=read.table(paste(fpx, "qc", sep="."), head=F)
    Ks=read.table(paste(fpx, "ks", sep="."), head=F)

    List=list()
    List1=list()
    List2=list()

    # Here analysis begins

    # RLE

    m=nrow(Rle)

    Mean=mean(Rle[,2])
    Sd=sd(Rle[,2])

    IQR_Mean=mean(Rle[,3])
    IQR_Sd=sd(Rle[,3])

    for (i in 1:m) {
        Name=as.character(Rle[i,1])
        # if (abs(Rle[i,2]-Mean) > (1.75 * Sd)) List[Name]=Name
        # if (abs(Rle[i,3]-IQR_Mean) > (1.75 * IQR_Sd)) List[Name]=Name
        if (abs(Rle[i,2]) > 0.075) List[Name]=Name
        if (abs(Rle[i,3]) > 0.75) List[Name]=Name
        if (abs(Rle[i,2]) > 0.075) List1[Name]=Name
        if (abs(Rle[i,3]) > 0.75) List1[Name]=Name
    }

    # NUSE

    Mean=mean(Nuse[,2])
    Sd=sd(Nuse[,2])

    IQR_Mean=mean(Nuse[,3])
    IQR_Sd=sd(Nuse[,3])

    for (i in 1:m) {
        Name=as.character(Rle[i,1])
        # if (abs(Nuse[i,2]-Mean) > (1.75 * Sd)) List[Name]=Name
        # if (abs(Nuse[i,3]-IQR_Mean) > (1.75 * IQR_Sd)) List[Name]=Name
        if (abs(Nuse[i,2]-1) > 0.075) List[Name]=Name
        if (abs(Nuse[i,3]) > 0.75) List[Name]=Name
        if (abs(Nuse[i,2]-1) > 0.075) List1[Name]=Name
        if (abs(Nuse[i,3]) > 0.75) List1[Name]=Name
    }

    # QC BioB

    for (i in 1:m) {
        Name=as.character(Rle[i,1])
        if (as.character(Qc[i,6]) == "A") List[Name]=Name
        if (as.character(Qc[i,6]) == "A") List2[Name]=Name
    }

    # QC scale
    Qc[,3]=abs(log2(Qc[,3]))

    for (i in 1:m) {
        Name=as.character(Rle[i,1])
        if (Qc[i,3] > 3) List[Name]=Name
        if (Qc[i,3] > 3) List2[Name]=Name
    }

    # And here we write

    Lv = as.vector(List)
    write.table(Lv, paste(fpx, "sum", sep="."), col.names=F, row.names=F, quote=F, sep="\n")
    Lv1 = as.vector(List1)
    write.table(Lv1, paste(fpx, "rmasum", sep="."), col.names=F, row.names=F, quote=F, sep="\n")
    Lv2 = as.vector(List2)
    write.table(Lv2, paste(fpx, "mas5sum", sep="."), col.names=F, row.names=F, quote=F, sep="\n")
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 1){
    compute_summary(d=args[1])
}
