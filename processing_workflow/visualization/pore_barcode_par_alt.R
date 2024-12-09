set=21
readalignment   ="no"

makegraphs      ="no"

expectbychance  ="no"

verbose         ="yes"

if(readalignment=="yes"){
suffix=paste("_11_17_seq_", set, ".pdf")
print("Reading in Alignment...")
bar.align=read.table("output_all_no_head.psla", sep="\t")
print("done")

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

print("Subsetting only 13+nt (inclusive) alignments...")
bar.align=bar.align[-which((bar.align$match)<13),]
print("done")

print("Ordering dataset...")
bar.align=bar.align[order(bar.align$Q_name),]
print("done")

bad.par=0
if (bad.par==1){
start=(set*50000)+1
end=(set+1)*50000
start.1=min(which(bar.align$Q_name==unique(bar.align$Q_name)[start]))
end.1  =max(which(bar.align$Q_name==unique(bar.align$Q_name)[end]))
#if (dim(bar.align)[1]>49999){
#bar.align=bar.align[start.1:end.1,]}
}

if (makegraphs=="yes"){
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
}
}


bar.align[,14] = sapply(bar.align[,14], as.character)

gaps=bar.align$T_gap_count>=1
bar.align.orig=bar.align
bar.align=bar.align[!gaps,]

mask1.m=matrix(nrow=length(unique(bar.align$Q_name)), ncol=8)

i=1
length=as.numeric(length(unique(bar.align$Q_name)))
readinloop=1
print("Beginning big loop") 
search.min=1
search.max=100
j=1
set.mask=1
#lapply(unique(bar.align$Q_name)[1:10000], function(ea){
for (ea in sort(unique(bar.align$Q_name))){
    ea = as.character(ea)
	print(ea)
	print(search.min)
	print(search.max)	
	if (search.min==Inf | as.numeric(search.min)<1){
		search.min=2
		search.max=1007071}
	print(ea)
	print(search.min)
	print(search.max)

	search.min=as.numeric(search.min)
	search.max=as.numeric(search.max)
    set.mask=which(bar.align$Q_name[search.min:search.max]==ea)    #indeces of the current read we are working on
	print(set.mask)	
	if (length(set.mask)<0){
		search.min=2
		search.max=1007071
		set.mask=which(bar.align$Q_name[search.min:search.max]==ea)
}


	search.min=min(set.mask)-2500
	search.max=search.min+5000	
		


	print(ea)
	print(search.min)
	print(search.max)
    bar.align.working=bar.align[set.mask,]  #subset on the current read we are working on

    max.val=max(bar.align.working$match)
    count=bar.align.working$match==max.val
    over.18=bar.align.working$match>=18     #18 or more nucleotides match in alignment out of a total of 24 nucleotides in barcode
    under.18=(bar.align.working$match<18 & bar.align.working$match>12)  #less than 18 nt match
    Q.start.early=bar.align.working$Q_start<101 | ((bar.align.working$Q_size - bar.align.working$Q_start) < 101 & (bar.align.working$Q_size - bar.align.working$Q_start) > 25)  #mask of indeces of alignments which start within 100 bp of read ends
    Q.start.late=bar.align.working$Q_start>100  | (bar.align.working$Q_size - bar.align.working$Q_start) > 100      #mask of indeces of alignments which do not start within 100 bp of read ends
    
    mask1=intersect(which(sapply(over.18, isTRUE)), which(sapply(Q.start.early, isTRUE)))   
    
    mask2=intersect(which(sapply(over.18, isTRUE)), which(sapply(Q.start.late, isTRUE)))
    
    mask3=intersect(which(sapply(under.18, isTRUE)), which(sapply(Q.start.early, isTRUE)))
    
    mask4=intersect(which(sapply(under.18, isTRUE)), which(sapply(Q.start.late, isTRUE)))
    
    bar.c=table(bar.align.working[count, 14])
    mask0=as.numeric(which.max(bar.c)[1])
    if (verbose=="yes"){
    print(paste("Processing read", readinloop, "out of", length, "reads"), sep=" ")
    print(ea)
    readinloop=readinloop+1}
    
    if (length(mask1) >0){
        mask1.m[i,1]=as.character(bar.align.working[mask1, 10])[1]
        mask1.m[i,2]=bar.align.working[mask1, 14][1]
        mask1.m[i,3]=bar.align.working[mask1, 12][1]
        mask1.m[i,4]=bar.align.working[mask1, 13][1]
        mask1.m[i,5]=bar.align.working[mask1, 11][1]
        mask1.m[i,6]=bar.align.working[mask1, 1][1]
        mask1.m[i,7]="mask1"
    }
    
    if (length(mask2) >0 & length(mask1) <= 0){
        mask1.m[i,1]=as.character(bar.align.working[mask2, 10])[1]
        mask1.m[i,2]=bar.align.working[mask2, 14][1]
        mask1.m[i,3]=bar.align.working[mask2, 12][1]
        mask1.m[i,4]=bar.align.working[mask2, 13][1]
        mask1.m[i,5]=bar.align.working[mask2, 11][1]
        mask1.m[i,6]=bar.align.working[mask2, 1][1]
        mask1.m[i,7]="mask2"
    }
    
    if (length(mask3) >0 & length(mask1) <= 0 & length(mask2) <= 0){
        mask1.m[i,1]=as.character(bar.align.working[mask3, 10])[1]
        mask1.m[i,2]=bar.align.working[mask3, 14][1]
        mask1.m[i,3]=bar.align.working[mask3, 12][1]
        mask1.m[i,4]=bar.align.working[mask3, 13][1]
        mask1.m[i,5]=bar.align.working[mask3, 11][1]
        mask1.m[i,6]=bar.align.working[mask3, 1][1]
        mask1.m[i,7]="mask3"
    }
    
    if (length(mask4) >0 & length(mask1) <= 0 & length(mask2) <= 0 & length(mask3) <= 0){
        mask1.m[i,1]=as.character(bar.align.working[mask4, 10])[1]
        mask1.m[i,2]=bar.align.working[mask4, 14][1]
        mask1.m[i,3]=bar.align.working[mask4, 12][1]
        mask1.m[i,4]=bar.align.working[mask4, 13][1]
        mask1.m[i,5]=bar.align.working[mask4, 11][1]
        mask1.m[i,6]=bar.align.working[mask4, 1][1]
        mask1.m[i,7]="mask4"
    }
    
    if (length(mask4) <= 0 & length(mask1) <= 0 & length(mask2) <= 0 & length(mask3) <= 0){
        mask1.m[i,1]=as.character(bar.align.working[mask0, 10])[1]
        mask1.m[i,2]=bar.align.working[mask0, 14][1]
        mask1.m[i,3]=bar.align.working[mask0, 12][1]
        mask1.m[i,4]=bar.align.working[mask0, 13][1]
        mask1.m[i,5]=bar.align.working[mask0, 11][1]
        mask1.m[i,6]=bar.align.working[mask0, 1][1]
        mask1.m[i,7]="mask0"
    }
    
mask1.m[i,8]=paste(mask1.m[i,1], "_", mask1.m[i,2], "_", mask1.m[i,3], "_", mask1.m[i,4], "_", mask1.m[i,5], "_", mask1.m[i,6], "_", mask1.m[i,7], ".fasta", sep="")
i=i+1    
}
#)   

print("Done with big loop...writing file...")
mask1.m=as.data.frame(mask1.m)
colnames(mask1.m)[1]="Name"
colnames(mask1.m)[2]="Barcode"
colnames(mask1.m)[3]="Query_Match_Start"
colnames(mask1.m)[4]="Query_Match_End"
colnames(mask1.m)[5]="Query_Length"
colnames(mask1.m)[6]="Match_num_nucl"
colnames(mask1.m)[7]="Mask_Num"
colnames(mask1.m)[8]="Filename"
write.csv(mask1.m, file=paste("align_index", print(Sys.time()), ".pdf", sep="_" ))

print("done")

if(expectbychance=="yes"){          #calculation of the probabilty of an alignment of a certain length given the average read length and bacode length, can use this to justify to lower nt length matches. This very simplistic model does imply 
pdf(file="expected_by_chance.pdf")
curve(1 - (1 - (1/(4^x)))^(2100), from=0, to=25, xlab="Number of Nucleotides matching Barcode", ylab="Probablity by Chance", main="P(N-mer match) in a 2.1 kb sequence by chance")
dev.off()}
