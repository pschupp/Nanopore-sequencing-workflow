library(data.table)
temp=0
if(temp==1){
suffix="_11_17_seq.pdf"
print("Reading in Alignment")
bar.align=read.table("output_all_no_head.psla", sep="\t")
print("Finished reading in Alignment")
colnames(bar.align)[1]="match"
colnames(bar.align)[2]="mismatch"
colnames(bar.align)[3]="rep.match"
colnames(bar.align)[4]="Ns"
colnames(bar.align)[5]="Q_gap_count"
colnames(bar.align)[6]="Q_gap_bases"
colnames(bar.align)[7]="T_gap_count"
colnames(bar.align)[8]="T_gap_bases"
colnames(bar.align)[9]="strand"
colnames(bar.align)[10]="Q_name"
colnames(bar.align)[11]="Q_size"
colnames(bar.align)[12]="Q_start"
colnames(bar.align)[13]="Q_end"
colnames(bar.align)[14]="T_name"
colnames(bar.align)[15]="T_size"
colnames(bar.align)[16]="T_start"
colnames(bar.align)[17]="T_end"
colnames(bar.align)[18]="Block_count"
colnames(bar.align)[19]="Block_sizes"
colnames(bar.align)[20]="qStarts"
colnames(bar.align)[21]="tStarts"
colnames(bar.align)[22]="qSeq"
colnames(bar.align)[23]="tSeq"

mask=(as.list(bar.align$Q_start)<100)
align=bar.align[mask,]
unalign=bar.align[!mask,]
exclusive.mask=unalign$Q_name %in% align$Q_name
unalign=unalign[!exclusive.mask,]
print("Making Graphs")
pdf(file=paste("Barcode_start", suffix))
hist(bar.align[bar.align$Q_start<2000,12], breaks=1000, xaxt="n", main="Start point of Barcodes", xlab="Start point of Barcode (bp)", ylab="Frequency")
axis(1, at = seq(0, 2000, by = 50), las=2)
dev.off()

pdf(paste(file="Barcode_start_tight", suffix)) 
hist(bar.align[bar.align$Q_start<100,12], breaks=1000, xaxt="n", main="Start point of Barcodes", xlab="Start point of Barcode (bp)", ylab="Frequency")
axis(1, at = seq(0, 100, by = 5), las=2)
dev.off()

pdf(paste(file="Barcode_match", suffix))
hist(bar.align$match, main="Nucleotide Matches of Barcodes", xlab="Number of Nucleotides Matched", ylab="Frequency", xaxt="n")
axis(1, at = seq(0, 25, by =1), las=1)
dev.off()

bar.align[,14] = sapply(bar.align[,14], as.character)

gaps=bar.align$T_gap_count>=1
bar.align.orig=bar.align
bar.align=bar.align[!gaps,]
}
mask1.m=matrix(nrow=length(unique(bar.align$Q_name)), ncol=8)

i=1
length=as.numeric(length(unique(bar.align$Q_name)))
readinloop=1
print("Beginning big loop") 

lapply(unique(bar.align$Q_name)[1:10000], function(ea){
#for (ea in unique(bar.align$Q_name)){
    ea = as.character(ea)
    set.mask=which(bar.align$Q_name==ea)
    max.val=max(bar.align$match[set.mask])
    count=bar.align$match[set.mask]==max.val
    over.18=bar.align$match[set.mask]>=18
    under.18=(bar.align$match[set.mask]<18 & bar.align$match[set.mask]>12)
    Q.start.early=bar.align$Q_start[set.mask]<101
    Q.start.late=bar.align$Q_start[set.mask]>100
    
    mask1=intersect(which(sapply(over.18, isTRUE)), which(sapply(Q.start.early, isTRUE)))
    
    mask2=intersect(which(sapply(over.18, isTRUE)), which(sapply(Q.start.late, isTRUE)))
    
    mask3=intersect(which(sapply(under.18, isTRUE)), which(sapply(Q.start.early, isTRUE)))
    
    mask4=intersect(which(sapply(under.18, isTRUE)), which(sapply(Q.start.late, isTRUE)))
    
    bar.c=table(bar.align[set.mask[count], 14])
    mask0=as.numeric(which.max(bar.c)[1])
    print(paste("Processing", readinloop, "out of", length), sep=" ")
    readinloop=readinloop+1
    print(ea)
    if (length(mask1) >0){
        mask1.m[i,1]=as.character(bar.align[set.mask[mask1], 10])[1]
        mask1.m[i,2]=bar.align[set.mask[mask1], 14][1]
        mask1.m[i,3]=bar.align[set.mask[mask1], 12][1]
        mask1.m[i,4]=bar.align[set.mask[mask1], 13][1]
        mask1.m[i,5]=bar.align[set.mask[mask1], 11][1]
        mask1.m[i,6]=bar.align[set.mask[mask1], 1][1]
        mask1.m[i,7]="mask1"
        i=i+1
    }
    
    if (length(mask2) >0 & length(mask1) <= 0){
        mask1.m[i,1]=as.character(bar.align[set.mask[mask2], 10])[1]
        mask1.m[i,2]=bar.align[set.mask[mask2], 14][1]
        mask1.m[i,3]=bar.align[set.mask[mask2], 12][1]
        mask1.m[i,4]=bar.align[set.mask[mask2], 13][1]
        mask1.m[i,5]=bar.align[set.mask[mask2], 11][1]
        mask1.m[i,6]=bar.align[set.mask[mask2], 1][1]
        mask1.m[i,7]="mask2"
        i=i+1
    }
    if (length(mask3) >0 & length(mask1) <= 0 & length(mask2) <= 0){
        mask1.m[i,1]=as.character(bar.align[set.mask[mask3], 10])[1]
        mask1.m[i,2]=bar.align[set.mask[mask3], 14][1]
        mask1.m[i,3]=bar.align[set.mask[mask3], 12][1]
        mask1.m[i,4]=bar.align[set.mask[mask3], 13][1]
        mask1.m[i,5]=bar.align[set.mask[mask3], 11][1]
        mask1.m[i,6]=bar.align[set.mask[mask3], 1][1]
        mask1.m[i,7]="mask3"
        i=i+1
    }
    if (length(mask4) >0 & length(mask1) <= 0 & length(mask2) <= 0 & length(mask3) <= 0){
        mask1.m[i,1]=as.character(bar.align[set.mask[mask4], 10])[1]
        mask1.m[i,2]=bar.align[set.mask[mask4], 14][1]
        mask1.m[i,3]=bar.align[set.mask[mask4], 12][1]
        mask1.m[i,4]=bar.align[set.mask[mask4], 13][1]
        mask1.m[i,5]=bar.align[set.mask[mask4], 11][1]
        mask1.m[i,6]=bar.align[set.mask[mask4], 1][1]
        mask1.m[i,7]="mask4"
        i=i+1
    }
    if (length(mask4) <= 0 & length(mask1) <= 0 & length(mask2) <= 0 & length(mask3) <= 0){
        mask1.m[i,1]=as.character(bar.align[set.mask[mask0], 10])[1]
        mask1.m[i,2]=bar.align[set.mask[mask0], 14][1]
        mask1.m[i,3]=bar.align[set.mask[mask0], 12][1]
        mask1.m[i,4]=bar.align[set.mask[mask0], 13][1]
        mask1.m[i,5]=bar.align[set.mask[mask0], 11][1]
        mask1.m[i,6]=bar.align[set.mask[mask0], 1][1]
        mask1.m[i,7]="mask0"
        i=i+1
    }
}
)   
for (j in c(1:dim(mask1.m)[1])){
    mask1.m[j,8]=paste(mask1.m[j,1], "_", mask1.m[j,2], "_", mask1.m[j,3], "_", mask1.m[j,4], "_", mask1.m[j,5], "_", mask1.m[j,6], "_", mask1.m[j,7], ".fasta", sep="")
} 
mask1.m=as.data.frame(mask1.m)
colnames(mask1.m)[1]="Name"
colnames(mask1.m)[2]="Barcode"
colnames(mask1.m)[3]="Query_Match_Start"
colnames(mask1.m)[4]="Query_Match_End"
colnames(mask1.m)[5]="Query_Length"
colnames(mask1.m)[6]="Match_num_nucl"
colnames(mask1.m)[7]="Mask_Num"
colnames(mask1.m)[8]="Filename"
write.csv(mask1.m, file="align.index.csv")

pdf(file="expected_by_chance.pdf")
expect=matrix(nrow=25, ncol=2)
for (k in c(1:25)){
expect[k,1]=k
expect[k,2]=1-(1-(1/(4^k)))^(2.1*1000)
}
plot(expect, xlab="Number of Nucleotides matching Barcode", ylab="Probablity by Chance", main="P(N-mer match) in a 2.1 kb sequence by chance")
dev.off()
