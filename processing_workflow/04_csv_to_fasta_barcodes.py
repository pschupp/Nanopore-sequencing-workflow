path="/Users/Patrick/Downloads/" # don't forget final "/"
file="cas_guides.csv"
wpath=paste(path, file, sep="")
barcodes=read.table(wpath, sep=",")
b.mat=matrix(dim(barcodes)[1]*2, ncol=1)

even=seq(2, dim(barcodes)[1]*2, by=2)
odd=seq(1,dim(barcodes)[1]*2, by=2)

names=barcodes[,1]
names=unlist(lapply(names, function(x)  paste(">", x, sep="")))
b.mat[odd]=names

seq=as.character(barcodes[,2])
b.mat[even]=seq

write.table(b.mat, file=paste(path, "singe_cell_barcodes.fasta", sep=""), row.names=F, col.names=F, quote=FALSE)