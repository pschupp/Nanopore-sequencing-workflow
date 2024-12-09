python -m HTSeq.scripts.count [options]  union --nonunique all -t exon <alignment_files> <gff_file>

#can also try in future more permissive -t gene so it can match anywhere.


/usr/local/bin/htseq-count -m union --nonunique all -a 0 -i gene_id --additional-attr=gene_name -t gene -s no ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/alignments_minimap2_o_min.sam ~/Documents/hg_align_db/GRCh37.p13_Gencode/gencode.v19.annotation.gff > htseq_counts.test

/usr/local/bin/htseq-count -m union --nonunique all -t exon -a 0 -i gene_id --additional-attr=gene_name -s no ~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/alignments_minimap2_o_min.sam ~/Documents/hg_align_db/GRCh37.p13_Gencode/gencode.v19.annotation.gff > htseq_counts_exon.test

/usr/local/bin/htseq-count -m union --nonunique all -t exon -a 0 -i gene_id --additional-attr=gene_name -s no ~/Documents/24_sec_seq/root/08_read_assigned_to_gene/split_BC25_test.sam ~/Documents/hg_align_db/GRCh37.p13_Gencode/gencode.v19.annotation.gff > htseq_counts_exon_split.test

library(data.table)
alig=fread(file="~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/alignments_minimap2_o_min.sam", sep="\t", skip=298, fill=TRUE, na.strings= "")
sub =grep('39', alig$V1)
header=read.table(file="/home/patrick/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/alignments_minimap2_o_min.sam", nrows=298 ,  sep="\n", header=FALSE)
alig=rbind(header, alig, fill=TRUE)

write.table(as.data.frame(alig)[sub,], file="split_BC25_test.sam", sep="\t", quote=FALSE, na="", row.names=FALSE, col.names=FALSE)

num=seq(25,49,1)
search=paste("_BC", num, sep="")

library(data.table)
alig=read.table(file="BC_all_reads_output_k14_ercc.sam", sep="\t", fill=TRUE, skip=298, na.strings= "",header=FALSE, col.names=1:23, stringsAsFactors=FALSE)
header=read.table(file="BC_all_reads_output_k14_ercc.sam", nrows=298 ,  sep="\n", header=FALSE, fill=TRUE, col.names=1:23, stringsAsFactors=FALSE)

for (ea in search){
sub=grep(ea, alig$X1)
alig.write=rbind(header, alig[sub,])
filename=paste("split", ea, "_alignment.sam", sep="")
write.table(as.data.frame(alig.write), file=filename, sep="\t", quote=FALSE, na="", row.names=FALSE, col.names=FALSE)
}

num=seq(25,49,1)
search=paste("_BC", num, sep="")

library(data.table)
alig=read.table(file="BC_all_reads_output_k14_38.p10.sam" , sep="\t", fill=TRUE, skip=556, na.strings= "",header=FALSE, col.names=1:23, stringsAsFactors=FALSE)
header=read.table(file="BC_all_reads_output_k14_38.p10.sam", nrows=556 ,  sep="\n", header=FALSE, fill=TRUE, col.names=1:23, stringsAsFactors=FALSE)

num=seq(25,49,1)
search=paste("_BC", num, sep="")
ercc.counts=as.data.frame(matrix(nrow=24, ncol=25))
ercc.counts[,1:2]=as.data.frame(table(alig$X3))
ercc.counts=ercc.counts[-1,]
i=3

for (ea in search){
sub=grep(ea, alig$X1)
alig.write=rbind(header, alig[sub,])
ercc.ind=which(alig.write$X2!=4)
alig.ercc=alig.write[ercc.ind,]
alig.ercc.counts=as.data.frame(table(alig.ercc$X3))

alig.ind=which(ercc.counts[,1] %in% alig.ercc.counts[,1])
ercc.counts[alig.ind,i]=try(alig.ercc.counts[,2])
i=i+1
}
ercc.counts[is.na(ercc.counts)]=0
rowSums(apply(as.matrix(ercc.counts[, 3:27]), c(1,2),as.numeric)) == ercc.counts[,2]
colnames(ercc.counts)[2]="rowSums"
colnames(ercc.counts)[3:ncol(ercc.counts)]=search
write.table(ercc.counts, file="ercc.counts.tsv", quote=FALSE, sep="\t", row.names=FALSE)


filename=paste("split", ea, "_alignment.sam", sep="")
write.table(as.data.frame(alig.write), file=filename, sep="\t", quote=FALSE, na="", row.names=FALSE, col.names=FALSE)



#version no ercc
library(data.table)
alig=read.table(file="BC_all_reads_output_k14_38.p10.primary.sam" , sep="\t", fill=TRUE, skip=556, na.strings= "",header=FALSE, col.names=1:23, stringsAsFactors=FALSE)
header=read.table(file="BC_all_reads_output_k14_38.p10.primary.sam", nrows=556 ,  sep="\n", header=FALSE, fill=TRUE, col.names=1:23, stringsAsFactors=FALSE)


num=seq(25,49,1)
search=paste("_BC", num, sep="")

i=3

for (ea in search){
sub=grep(ea, alig$X1)
alig.write=rbind(header, alig[sub,])

i=i+1

filename=paste("split", ea, "_alignment.sam", sep="")
write.table(as.data.frame(alig.write), file=filename, sep="\t", quote=FALSE, na="", row.names=FALSE, col.names=FALSE)
}

#version no ercc, no header
library(data.table)
alig=read.table(file="BC_all_reads_output_k14_38.p10.primary.sam" , sep="\t", fill=TRUE, na.strings= "",header=FALSE, col.names=1:23, stringsAsFactors=FALSE)



num=seq(25,49,1)
search=paste("_BC", num, sep="")

i=3

for (ea in search){
sub=grep(ea, alig$X1)
alig.write=alig[sub,]

i=i+1

filename=paste("split", ea, "_alignment.sam", sep="")
write.table(as.data.frame(alig.write), file=filename, sep="\t", quote=FALSE, na="", row.names=FALSE, col.names=FALSE)
}