search=list.files(pattern="alignment.sam.counts.txt")

input.s=read.table(search[1], fill=TRUE, stringsAsFactors=FALSE)
input.s=input.s[1:(dim(input.s)[1]-5), ]    #last 5 rows have statistics on non-alignments: no feature, ambigious, not aligned, etc.

for (ea in search[2:length(search)]){
input=read.table(ea, fill=TRUE, stringsAsFactors=FALSE)
input=input[1:(dim(input)[1]-5), ]
input.s=cbind(input.s, input[,3])
}

colnames(input.s)[1:2]=c("gene_id", "gene_name")
colnames(input.s)[3:length(colnames(input.s))]=sapply(search, substr, 7, 10)

colsums=apply(input.s[,3:dim(input.s)[2]],2, sum)
rowsums=apply(input.s[,3:dim(input.s)[2]],1, sum)
sum(colsums)

unique =as.data.frame(table(input.s$gene_name))
unique.gene=unique[which(unique[,2]>1),1]   #which gene names occur more than once
index.rm=c()
for (gene in unique.gene){
	index.w=which(input.s$gene_name==gene)
	input.w=input.s[index.w, ]
	colsums.w=apply(input.w[,3:ncol(input.w)], 2, sum)
	union.w=c(input.s[index.w[1], 1:2], colsums.w)
	input.s=rbind(input.s, union.w)
	index.rm=c(index.rm, index.w)
}
input.s=input.s[-index.rm,]
rowsums=apply(input.s[,3:dim(input.s)[2]],1, sum)
input.s=input.s[-which(rowsums==0),]

colsums=apply(input.s[,3:dim(input.s)[2]],2, sum)
sum(colsums)

write.table(input.s, file="nano_24_sec_GRCh38_alig_no_ercc.tsv", quote=FALSE, sep="\t", row.names=FALSE)

ercc.counts=read.table("~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh37.p7_alignment/ERCC_alignment/ercc.counts.tsv", sep="\t")
ercc.counts=ercc.counts[,-c(10,18)] #remove BC32 and 40 which were not used
colnames(ercc.counts)=colnames(input.s)
ercc.counts[,2]=ercc.counts[,1]
input.s.ercc=rbind(input.s, ercc.counts)

write.table(input.s.ercc, file="nano_24_sec_GRCh38_alig_ercc.tsv", quote=FALSE, sep="\t", row.names=FALSE)
rownames(input.s.ercc)=input.s.ercc$gene_name
input.s.ercc=input.s.ercc[,-c(1)]
colnames(input.s.ercc)[1]="Gene"
zeros=apply(input.s.ercc[,c(2:ncol(input.s.ercc))], 1, function(x) sum(x==0))
table(zeros)

write.table(input.s.ercc, file="nano_24_sec_working.tsv", quote=FALSE, sep="\t", row.names=FALSE)
many_zeros=which(zeros>17)
input.s.ercc.min=input.s.ercc[-many_zeros,]
write.table(input.s.ercc.min, file="nano_24_sec_min.tsv", quote=FALSE, sep="\t", row.names=FALSE)

expr=read.table("nano_24_sec_min.tsv", stringsAsFactors=FALSE, sep="\t", header=TRUE)
expr=expr[-5683,]
expr.n=t(sapply(expr[,c(2:ncol(expr))],as.numeric))
colnames(expr.n)=expr[,1]
nano.pca=prcomp(expr.n, scale.=TRUE, center=TRUE)
nano.pca.eig=nano.pca$sdev^2
nano.pca.pvar=nano.pca.eig/sum(nano.pca.eig)
xlab.char=paste("PC1 ", signif(nano.pca.pvar[1]*100, digits=3), " % var. explained", sep="")
ylab.char=paste("PC2 ", signif(nano.pca.pvar[2]*100, digits=3), " % var. explained", sep="")

pdf(file="nanopore_24s_expression_norm_log2_pca.pdf")
nano.pca.plot=biplot(nano.pca, ylabs=rep("", nrow(expr)), var.axes=FALSE, main="Plot of First 2 PCs of Log2-Norm Data", xlab=xlab.char, ylab=ylab.char )
dev.off()

pdf(file="nanopore_24s_expression_norm_log2_pca_alt.pdf")
nano.pca.plot=biplot(nano.pca, ylabs=rep(".", nrow(expr)), var.axes=TRUE, main="Plot of First 2 PCs of Log2-Norm Data", xlab=xlab.char, ylab=ylab.char )
dev.off()

pdf(file="nanopore_24s_expression_norm_log2_pca_alt2.pdf")
nano.pca.plot=biplot(nano.pca, ylabs=rep(".", nrow(expr)), var.axes=FALSE, main="Plot of First 2 PCs of Log2-Norm Data", xlab=xlab.char, ylab=ylab.char )
dev.off()

expr=read.table("nano_24_sec_GRCh38_alig_no_ercc.tsv", stringsAsFactors=FALSE, sep="\t", header=TRUE)
expr=expr[-1]
zeros=apply(expr[,c(2:ncol(expr))], 1, function(x) sum(x==0))
table(zeros)
many_zeros=which(zeros>17)
expr=expr[-many_zeros,]
rownames(expr)=expr[,1]
expr=expr[,-1]


expr.n=t(expr[,-1])


nano.pca=prcomp(expr.n, scale.=TRUE, center=TRUE)
nano.pca.eig=nano.pca$sdev^2
nano.pca.pvar=nano.pca.eig/sum(nano.pca.eig)
xlab.char=paste("PC1 ", signif(nano.pca.pvar[1]*100, digits=3), " % var. explained", sep="")
ylab.char=paste("PC2 ", signif(nano.pca.pvar[2]*100, digits=3), " % var. explained", sep="")

pdf(file="nanopore_24s_expression_norm_log2_pca.pdf")
nano.pca.plot=biplot(nano.pca, ylabs=rep("", nrow(expr)), var.axes=FALSE, main="Plot of First 2 PCs of Log2-Norm Data", xlab=xlab.char, ylab=ylab.char )
dev.off()

pdf(file="nanopore_24s_expression_norm_log2_pca_alt.pdf")
nano.pca.plot=biplot(nano.pca, ylabs=rep(".", nrow(expr)), var.axes=TRUE, main="Plot of First 2 PCs of Log2-Norm Data", xlab=xlab.char, ylab=ylab.char )
dev.off()

pdf(file="nanopore_24s_expression_norm_log2_pca_alt2.pdf")
nano.pca.plot=biplot(nano.pca, ylabs=rep(".", nrow(expr)), var.axes=FALSE, main="Plot of First 2 PCs of Log2-Norm Data", xlab=xlab.char, ylab=ylab.char )
dev.off()


nano.ruvg2.pca=prcomp(t(nano.ruvg2.df), center=TRUE, scale.=TRUE)
nano.ruvg2.pca.eig=nano.ruvg2.pca$sdev^2
nano.ruvg2.pca.pvar=nano.ruvg2.pca.eig/sum(nano.ruvg2.pca.eig)
xlab.char=paste("PC1 ", signif(nano.ruvg2.pca.pvar[1]*100, digits=3), " % var. explained", sep="")
ylab.char=paste("PC2 ", signif(nano.ruvg2.pca.pvar[2]*100, digits=3), " % var. explained", sep="")

pdf(file=paste("nanopore_24s_expression_norm_log2_ruvg_pca_k", k.ruv, ".pdf", sep=""))
nano.pca.plot=biplot(nano.ruvg2.pca, ylabs=rep("", nrow(expr.norm)), var.axes=FALSE, main="Plot of First 2 PCs of Log2-Norm Data", xlab=xlab.char, ylab=ylab.char )
dev.off()

pdf(file=paste("nanopore_24s_expression_norm_log2_ruvg_pca_alt_k", k.ruv, ".pdf", sep=""))
nano.pca.plot=biplot(nano.ruvg2.pca, ylabs=rep(".", nrow(expr.norm)), var.axes=TRUE, main="Plot of First 2 PCs of Log2-Norm Data", xlab=xlab.char, ylab=ylab.char )
dev.off()

pdf(file=paste("nanopore_24s_expression_norm_log2_ruvg_pca_alt2_k", k.ruv, ".pdf", sep=""))
nano.pca.plot=biplot(nano.ruvg2.pca, ylabs=rep(".", nrow(expr.norm)), var.axes=FALSE, main="Plot of First 2 PCs of Log2-Norm Data", xlab=xlab.char, ylab=ylab.char )
dev.off()
}


#gonna also try all the above with log2(x+1) data, just to make sure that it results in the same values

####don't need this, ruvg uses raw count data anyway#####
#input.s.norm=input.s
#for (col in seq(3, dim(input.s)[2])){
#	input.s.norm[,col]=input.s[,col]/colsums[col-2]
#}
#input.s.norm=input.s
#for (col in seq(3, dim(input.s)[2])){
#	input.s.norm[,col]=input.s[,col]/colsums[col-2]
#}