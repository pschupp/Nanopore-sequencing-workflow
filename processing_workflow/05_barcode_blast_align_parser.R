rm(list=ls())
require(data.table)
path="/home/patrick/Documents/24_sec_seq/root/03_basecalled_read_fasta/"
file="barcode_output.pslx"
print("Reading in Files")
input=paste(path, file, sep="")
bar.align=fread(input, sep='\t')

print("Finished Reading in Files")
colnames(bar.align)[1]="match"
colnames(bar.align)[2]="mismatch"
colnames(bar.align)[3]="Q_seq"
colnames(bar.align)[4]="T_seq"
colnames(bar.align)[5]="Q_gap_count"
colnames(bar.align)[6]="Q_gap_total"
colnames(bar.align)[7]="ppos"
colnames(bar.align)[8]="qframe"
colnames(bar.align)[9]="strand"
colnames(bar.align)[10]="Q_name"
colnames(bar.align)[11]="Q_size"
colnames(bar.align)[12]="Q_start"
colnames(bar.align)[13]="Q_end"
colnames(bar.align)[14]="T_name"
colnames(bar.align)[15]="T_size"
colnames(bar.align)[16]="T_start"
colnames(bar.align)[17]="T_end"
colnames(bar.align)[18]="Length Alignment"
colnames(bar.align)[19]="Score"
colnames(bar.align)[20]="bitscore"
colnames(bar.align)[21]="sframe"
colnames(bar.align)[22]="frames"
colnames(bar.align)[23]="evalue"

bar.align=bar.align[order(bar.align$Q_name),]

matches.end=matrix(0, ncol=1, nrow=dim(bar.align)[1])

print('Finished Reading and Processing Input Files')

mask1.m=matrix(nrow=length(unique(bar.align$Q_name)), ncol=7)
doubles.m=matrix(ncol=23)
names(doubles.m)=names(bar.align)
i=1
length=as.numeric(length(unique(bar.align$Q_name)))
readinloop=1
print("Beginning big loop") 
search.min=1
search.max=10000
j=1
set.mask=1
#lapply(unique(bar.align$Q_name)[1:10000], function(ea){

pb <- txtProgressBar(min = 1, max = length(sort(unique(bar.align$Q_name))), style = 3)

for (ea in sort(unique(bar.align$Q_name))){
    ea = as.character(ea)
    write=1
    set.mask=which(bar.align$Q_name[search.min:search.max]==ea)+search.min-1   #indeces of the current read we are working on
#
    bar.align.working=as.data.frame(bar.align[set.mask,])  #subset on the current read we are working on
	search.min=set.mask[length(set.mask)]
	search.max=search.min+10000
    bar.align.working.match=bar.align.working$match

    max_match=max(bar.align.working.match)
    which_max=which(bar.align.working.match==max_match)
    ends=which(bar.align.working$Q_start %in% c(seq(1,100,1), seq(bar.align.working$Q_size[1]-100,bar.align.working$Q_size[1],1)))
    index.win=intersect(which_max, ends)
    if (length(unique(bar.align.working$T_name[index.win]))< length(bar.align.working$T_name[index.win])){
        bar.win=rownames(data.frame(table(bar.align.working$T_name[index.win])>1))[1]
        index.win=which(bar.align.working$T_name==bar.win)[1]
    }
    if (length(index.win)<1){
        avg_start=mean(bar.align.working[order(bar.align.working$match, decreasing=TRUE)[1:1^ceiling(log10((dim(bar.align.working)[1]/10)))],12])
        win2=which.min(abs(bar.align.working$Q_start[which_max]-avg_start))
        index.win=which_max[win2]
    }
    if(length(index.win)>1){
        if (sum(bar.align.working$evalue[index.win]<0.1)>1){
            index.win=which(bar.align.working$evalue[index.win]<0.1)[1]    
        } else {
        doubles.m=rbind(doubles.m, as.matrix(bar.align.working))
        write=0
        }
        
    }    
    if (write==1){
        mask1.m[i,1]=as.character(bar.align.working[index.win, 10])[1]
        mask1.m[i,2]=as.character(bar.align.working[index.win, 14][1])
        mask1.m[i,3]=as.numeric(bar.align.working[index.win, 12][1])
        mask1.m[i,4]=as.numeric(bar.align.working[index.win, 13][1])
        mask1.m[i,5]=as.numeric(bar.align.working[index.win, 11][1])
        mask1.m[i,6]=as.numeric(bar.align.working[index.win, 1][1])
        mask1.m[i,7]=paste(mask1.m[i,1], "_", mask1.m[i,2], "_", mask1.m[i,3], "_", mask1.m[i,4], "_", mask1.m[i,5], "_", mask1.m[i,6],  ".fasta", sep="")
    }
    setTxtProgressBar(pb, i)
    i=i+1    
}   
close(pb)
print("Done with big loop...writing file...")
mask1.m=as.data.frame(mask1.m)
colnames(mask1.m)[1]="Name"
colnames(mask1.m)[2]="Barcode"
colnames(mask1.m)[3]="Query_Match_Start"
colnames(mask1.m)[4]="Query_Match_End"
colnames(mask1.m)[5]="Query_Length"
colnames(mask1.m)[6]="Match_num_nucl"
colnames(mask1.m)[7]="Filename"

write.csv(mask1.m, file=paste(path, "barcode_align_analysis.csv", sep=""))
doubles.m=as.data.frame(doubles.m)
write.csv(doubles.m, file=paste(path, "barcode_align_doubles.csv", sep=""))

print("done")