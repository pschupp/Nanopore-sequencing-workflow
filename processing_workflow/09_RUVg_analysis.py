library('RUVSeq')
library('WGCNA'

###f(x) input
project_name="nano_24_sec_seq"
expr.counts=read.table("nano_24_sec_min.tsv",  sep="\t", header=TRUE, stringsAsFactors=FALSE)
expr.counts=expr.counts[-seq(5683, 5690),]
ercc.counts=read.table("~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh37.p7_alignment/ERCC_alignment/ercc.counts.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
oldham.genes=read.table("Oldham_fidelity_4_celltypes_genes.csv",  sep=",", header=TRUE, stringsAsFactors=FALSE)
oldham.genes=oldham.genes[,-1]
#processing with only numbers in actual dataframe for easier numerical processing later on. rownames are genes/ercc, colnames are samples
rownames(expr.counts)=expr.counts[,1]
expr.counts=expr.counts[,-1]
genes=rownames(expr.counts)
expr.counts=as.data.frame(lapply(expr.counts, as.numeric))
rownames(expr.counts)=genes

expr.colsums=apply(expr.counts,2, sum)
expr.norm=as.data.frame(matrix(nrow=nrow(expr.counts), ncol=ncol(expr.counts)))
for(col in seq(1,ncol(expr.counts))){
    expr.norm[,col]=expr.counts[,col]/expr.colsums[col]
}
expr.norm.l2=log2(expr.norm+1)
colnames(expr.norm.l2)=colnames(expr.counts)
rownames(expr.norm.l2)=rownames(expr.counts)

rownames(ercc.counts)=ercc.counts[,1]
ercc.counts=ercc.counts[,-c(1,2)]
ercc.counts=ercc.counts[,-c(8,16)]
colnames(ercc.counts)=colnames(expr.counts)
ercc.counts=ercc.counts[-which(apply(ercc.counts, 1, function(x) sum(x==0))>12),] #remove ercc synthetic transcripts for which there is a darth of count data

#using ruvg we are measuring two metrics to validate our approach. first we will be using PCA analysis to look at how clustering improves through RUVg normaliztion with the aim to delineate clear outliers, if applicable, and reduce the variation between samples dervied from non-biological sources. second we will use the log2 ratio of gene count / median gene count over all sections. RUVg should standardize the distribution of these values across all sections.

#first we will plot these plots with baseline data, log2 normalized and not log2 normalized.

expr.ruvg=RUVg(x=as.matrix(expr.counts), cIdx=as.matrix(ercc.counts), k=1, isLog=FALSE)
expr.ruvg.df=as.data.frame(expr.ruvg$normalizedCounts)
expr.ruvg.pca=prcomp(t(expr.ruvg.df), center=TRUE, scale.=TRUE)
expr.ruvg.pca.eig=expr.ruvg.pca$sdev^2
expr.ruvg.pca.pvar=expr.ruvg.pca.eig/sum(expr.ruvg.pca.eig)
xlab.char=paste("PC1 ", signif(expr.ruvg.pca.pvar[1]*100, digits=3), " % var. explained", sep="")
ylab.char=paste("PC2 ", signif(expr.ruvg.pca.pvar[2]*100, digits=3), " % var. explained", sep="")

pdf(file=paste(project_name, "pca_count_k0.pdf", sep="_"))
expr.pca.plot=biplot(expr.ruvg.pca, ylabs=rep("", nrow(expr.counts)), var.axes=FALSE, main="Plot of First 2 PCs of Log2-Norm Data", xlab=xlab.char, ylab=ylab.char )
dev.off()

pdf(file=paste(project_name, "pca_vec_count_k0.pdf", sep="_"))
expr.pca.plot=biplot(expr.ruvg.pca, ylabs=rep(".", nrow(expr.counts)), var.axes=FALSE, main="Plot of First 2 PCs of Log2-Norm Data", xlab=xlab.char, ylab=ylab.char )
dev.off()

expr.rmedian=as.numeric(apply((expr.counts+1),1, median))
expr.counts.r=as.data.frame(matrix(nrow=nrow(expr.counts),ncol=ncol(expr.counts)))
colnames(expr.counts.r)=colnames(expr.counts)
rownames(expr.counts.r)=rownames(expr.counts)

for (row in seq(1,nrow(expr.counts))){
    expr.counts.r[row,]=lapply(expr.counts[row,]+1, function(x) x/expr.rmedian[row])
}
expr.counts.r.l2=log2(expr.counts.r)
expr.counts.r.0=expr.counts.r.l2
expr.counts.r[expr.counts.r==0]=NA
expr.counts.r.l2[expr.counts.r.l2==0]=NA


pdf(file=paste(project_name, "rle_count_k0.pdf", sep="_"), width=18)
boxplot(expr.counts.r.l2, main="RLE of Unnormalized Data", ylab="RLE log2(count/median count)", xlab="Samples", notch=TRUE)
dev.off()

pdf(file=paste(project_name, "re_count_k0.pdf", sep="_"), width=18)
boxplot(expr.counts.r, main="RE of Unnormalized Data", ylab="RE (count/median count)", xlab="Samples", notch=TRUE)
dev.off()

pdf(file=paste(project_name, "rle_count_inc0_k0.pdf", sep="_"), width=18)
boxplot(expr.counts.r.0, main="RLE of Unnormalized Data", ylab="RLE log2(count/median count)", xlab="Samples", notch=TRUE)
dev.off()

k.ruv=0
write.table(expr.counts, quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t", file=paste(project_name, "_expression_matrix_log2_norm_k", k.ruv, ".tsv", sep=""))

###now doing iterative RUVg with various k's

for (k.ruv in c(seq(1, nrow(ercc.counts)))){
expr.ruvg=RUVg(x=as.matrix(expr.counts), cIdx=as.matrix(ercc.counts), k=k.ruv, isLog=FALSE)
expr.ruvg.df=as.data.frame(expr.ruvg$normalizedCounts)
expr.ruvg.df=log2(expr.ruvg.df+1)

write.table(expr.ruvg.df, quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t", file=paste(project_name, "_expression_matrix_log2_norm_k", k.ruv, ".tsv", sep=""))

vars=apply(expr.ruvg.df, 1, var)
expr.ruvg.df=expr.ruvg.df[which(vars!=0),]
if (dim(expr.ruvg.df)[1]>0){
expr.ruvg.pca=prcomp(t(expr.ruvg.df), center=TRUE, scale.=TRUE)
expr.ruvg.pca.eig=expr.ruvg.pca$sdev^2
expr.ruvg.pca.pvar=expr.ruvg.pca.eig/sum(expr.ruvg.pca.eig)
xlab.char=paste("PC1 ", signif(expr.ruvg.pca.pvar[1]*100, digits=3), " % var. explained", sep="")
ylab.char=paste("PC2 ", signif(expr.ruvg.pca.pvar[2]*100, digits=3), " % var. explained", sep="")

pdf(file=paste(project_name, "_pca_count_k", k.ruv, ".pdf", sep=""))
expr.pca.plot=biplot(expr.ruvg.pca, ylabs=rep("", nrow(expr.ruvg.df)), var.axes=FALSE, main="Plot of First 2 PCs of Log2-Norm Data", xlab=xlab.char, ylab=ylab.char )
dev.off()

pdf(file=paste(project_name, "_pca_vec_count_k", k.ruv, ".pdf", sep=""))
expr.pca.plot=biplot(expr.ruvg.pca, ylabs=rep(".", nrow(expr.ruvg.df)), var.axes=FALSE, main="Plot of First 2 PCs of Log2-Norm Data", xlab=xlab.char, ylab=ylab.char )
dev.off()

expr.rmedian=as.numeric(apply((expr.ruvg.df+1),1, median))
expr.counts.r=as.data.frame(matrix(nrow=nrow(expr.ruvg.df),ncol=ncol(expr.ruvg.df)))
colnames(expr.counts.r)=colnames(expr.ruvg.df)
rownames(expr.counts.r)=rownames(expr.ruvg.df)

for (row in seq(1,nrow(expr.ruvg.df))){
    expr.counts.r[row,]=lapply(expr.ruvg.df[row,]+1, function(x) x/expr.rmedian[row])
}
expr.counts.r.l2=log2(expr.counts.r)
expr.counts.r.l2[expr.counts.r.l2==0]=NA


pdf(file=paste(project_name, "_rle_count_k", k.ruv, ".pdf", sep=""), width=18)
boxplot(expr.counts.r.l2, main="RLE of Unnormalized Data", ylab="RLE log2(count/median count)", xlab="Samples", notch=TRUE)
dev.off()
}
}

##need to work automating this
k0=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k0.tsv")
k1=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k1.tsv")
k2=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k2.tsv")
k3=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k3.tsv")
k4=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k4.tsv")
k5=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k5.tsv")
k6=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k6.tsv")
k7=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k7.tsv")
k8=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k8.tsv")
k9=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k9.tsv")
k10=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k10.tsv")
k11=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k11.tsv")
k12=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k12.tsv")
k13=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k13.tsv")
k14=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k14.tsv")
k15=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k15.tsv")
k16=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k16.tsv")
k17=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k17.tsv")
k18=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k18.tsv")
k19=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k19.tsv")
k20=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k20.tsv")
k21=read.table(header=TRUE, sep="\t", file="nano_24_sec_seq_expression_matrix_log2_norm_k21.tsv")

expr.counts=log2(expr.counts+1)

for (type in seq(1,4,1)){
    expr.ind=which(rownames(expr.counts) %in% oldham.genes[,type])
    expr.norm.ind=which(rownames(expr.norm.l2) %in% oldham.genes[,type])
    k1.ind=which(rownames(k1) %in% oldham.genes[,type])
    k2.ind=which(rownames(k2) %in% oldham.genes[,type])
    k3.ind=which(rownames(k3) %in% oldham.genes[,type])
    k4.ind=which(rownames(k4) %in% oldham.genes[,type])
    k5.ind=which(rownames(k5) %in% oldham.genes[,type])
    k6.ind=which(rownames(k6) %in% oldham.genes[,type])
    k7.ind=which(rownames(k7) %in% oldham.genes[,type])
    k8.ind=which(rownames(k8) %in% oldham.genes[,type])
    k9.ind=which(rownames(k9) %in% oldham.genes[,type])
    k10.ind=which(rownames(k10) %in% oldham.genes[,type])
    k11.ind=which(rownames(k11) %in% oldham.genes[,type])
    k12.ind=which(rownames(k12) %in% oldham.genes[,type])
    k13.ind=which(rownames(k13) %in% oldham.genes[,type])
    k14.ind=which(rownames(k14) %in% oldham.genes[,type])
    k15.ind=which(rownames(k15) %in% oldham.genes[,type])
    k16.ind=which(rownames(k16) %in% oldham.genes[,type])
    k17.ind=which(rownames(k17) %in% oldham.genes[,type])
    k18.ind=which(rownames(k18) %in% oldham.genes[,type])
    k19.ind=which(rownames(k19) %in% oldham.genes[,type])
    k20.ind=which(rownames(k20) %in% oldham.genes[,type])
    k21.ind=which(rownames(k21) %in% oldham.genes[,type])
    
    expr.bicor=bicor(expr.counts[expr.ind,])
    expr.bicor=expr.bicor[upper.tri(expr.bicor)]
    expr.norm.bicor=bicor(expr.norm.l2[expr.norm.ind,])
    expr.norm.bicor=expr.norm.bicor[upper.tri(expr.norm.bicor)]
    
    k1.bicor=bicor(k1[k1.ind,])
    k1.bicor=k1.bicor[upper.tri(k1.bicor)]
    k2.bicor=bicor(k2[k2.ind,])
    k2.bicor=k2.bicor[upper.tri(k2.bicor)]
    k3.bicor=bicor(k3[k3.ind,])
    k3.bicor=k3.bicor[upper.tri(k3.bicor)]
    k4.bicor=bicor(k4[k4.ind,])
    k4.bicor=k4.bicor[upper.tri(k4.bicor)]
    k5.bicor=bicor(k5[k5.ind,])
    k5.bicor=k5.bicor[upper.tri(k5.bicor)]
    k6.bicor=bicor(k6[k6.ind,])
    k6.bicor=k6.bicor[upper.tri(k6.bicor)]
    k7.bicor=bicor(k7[k7.ind,])
    k7.bicor=k7.bicor[upper.tri(k7.bicor)]
    k8.bicor=bicor(k8[k8.ind,])
    k8.bicor=k8.bicor[upper.tri(k8.bicor)]
    k9.bicor=bicor(k9[k9.ind,])
    k9.bicor=k9.bicor[upper.tri(k9.bicor)]
    k10.bicor=bicor(k10[k10.ind,])
    k10.bicor=k10.bicor[upper.tri(k10.bicor)]
    k11.bicor=bicor(k11[k11.ind,])
    k11.bicor=k11.bicor[upper.tri(k11.bicor)]
    k12.bicor=bicor(k12[k12.ind,])
    k12.bicor=k12.bicor[upper.tri(k12.bicor)]
    k13.bicor=bicor(k13[k13.ind,])
    k13.bicor=k13.bicor[upper.tri(k13.bicor)]
    k14.bicor=bicor(k14[k14.ind,])
    k14.bicor=k14.bicor[upper.tri(k14.bicor)]
    k15.bicor=bicor(k15[k15.ind,])
    k15.bicor=k15.bicor[upper.tri(k15.bicor)]
    k16.bicor=bicor(k16[k16.ind,])
    k16.bicor=k16.bicor[upper.tri(k16.bicor)]
    k17.bicor=bicor(k17[k17.ind,])
    k17.bicor=k17.bicor[upper.tri(k17.bicor)]
    k18.bicor=bicor(k18[k18.ind,])
    k18.bicor=k18.bicor[upper.tri(k18.bicor)]
    k19.bicor=bicor(k19[k19.ind,])
    k19.bicor=k19.bicor[upper.tri(k19.bicor)]
    k20.bicor=bicor(k20[k20.ind,])
    k20.bicor=k20.bicor[upper.tri(k20.bicor)]
    k21.bicor=bicor(k21[k21.ind,])
    k21.bicor=k21.bicor[upper.tri(k21.bicor)]
    
    
    type.c=colnames(oldham.genes)[type]
    type.cor=data.frame(matrix(nrow=50, ncol=11))
    type.cor=cbind(expr.bicor, expr.norm.bicor, k1.bicor, k2.bicor, k3.bicor, k4.bicor, k5.bicor, k6.bicor, k7.bicor, k8.bicor, k9.bicor, k10.bicor, k11.bicor, k12.bicor, k13.bicor, k14.bicor, k15.bicor, k16.bicor, k17.bicor, k18.bicor, k19.bicor, k20.bicor, k21.bicor)
    colnames(type.cor)=c("Unnormalized", "Global-Scaling", "RUVg, k=1", "RUVg, k=2", "RUVg, k=3", "RUVg, k=4", "RUVg, k=5", "RUVg, k=6", "RUVg, k=7", "RUVg, k=8", "RUVg, k=9", "RUVg, k=10", "RUVg, k=11", "RUVg, k=12", "RUVg, k=13", "RUVg, k=14", "RUVg, k=15", "RUVg, k=16", "RUVg, k=17", "RUVg, k=18", "RUVg, k=19", "RUVg, k=20", "RUVg, k=21")
    pdf(file=paste("bwmcor_top50", type.c, "genes.pdf", sep="_"), width=15)
    boxplot(type.cor, ylab="Biweight Midcorrelation", main=paste("Biweight Midcor of Top 50", type.c, "genes", sep=" "))
    dev.off()
    
}