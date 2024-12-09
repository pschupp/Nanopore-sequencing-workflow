library(WVPlots)

meta=read.table("sequencing_summary.txt", sep="\t")
colnames(meta) <- sapply(meta[1,] ,as.character)
meta=meta[-1,]
meta=meta[-c(which(meta[,1]=="filename")),]

suffix="_11_17_seq.pdf"
suffix2="_sp_1800"


meta$sequence_length_template=as.numeric(as.character(meta$sequence_length_template))
meta$mean_qscore_template=as.numeric(as.character(meta$mean_qscore_template))
meta$duration=as.numeric(as.character(meta$duration))


length.mean=mean(meta$sequence_length_template)
length.sd=sd(meta$sequence_length_template)
length.mask=meta$sequence_length_template<(length.mean+6*length.sd)
if (is.na(table(length.mask)[2])==FALSE){
if (table(length.mask)[1]/table(length.mask)[2]>0.0002){
    warning("Unusually high fraction (>0.02%) of relatively long (>6*SD) reads", immediate=TRUE)
}}

pdf(file=paste("Read_Length_v_Quality_smoothscatter", suffix))
try(smoothScatter(meta$mean_qscore_template[length.mask], meta$sequence_length_template[length.mask], nbin=2000, xlab= "Read Q-Score", ylab="Read Length (bp)", main="Read Length and Q-Score"))
dev.off()

pdf(file=paste("Histogram_Sequence_Length", suffix))
hist(meta$sequence_length_template[length.mask], xlab="Read Length (bp)", ylab="Count", main="Distribution of Read Lengths")
dev.off()

pdf(file=paste("Histogram_Q_Scores", suffix))
hist(meta$mean_qscore_template[length.mask], xlab="Q-Score", ylab="Count", main="Distribution of Q-Scores")
dev.off()

length.mean.d=mean(meta$duration)
length.sd.d=sd(meta$duration)
length.mask.d=meta$duration<(length.mean.d+6*length.sd.d)
if (is.na(table(length.mask.d)[2])==FALSE){
    if (table(length.mask.d)[1]/table(length.mask.d)[2]>0.001){
    warning("Unusually high fraction (> 0.1%)of relatively high duration (>6*SD) reads", immediate=TRUE)
}}

pdf(file=paste("Histogram_Read_Duration", suffix))
hist(meta$duration[length.mask.d], xlab="Duration of Reads", ylab="Count", main="Distribution of Read Durations")
dev.off()

library(WVPlots)

meta.samp=meta[sample(c(1:dim(meta)[1]), 4000),]
if (4000/dim(meta)[1]<0.001){
    warning("Subsampling less than 0.1% of sample for Q_Score_Read_Length Plot which may not be representative", immediate=TRUE)
}

colnames(meta.samp)[14]="Q_Score"
colnames(meta.samp)[13]="Read_Length"

meta.samp[,14] <- sapply(meta.samp[,14] ,as.character)
meta.samp[,14] <- sapply(meta.samp[,14] ,as.numeric)

meta.samp[,13] <- sapply(meta.samp[,13] ,as.character)
meta.samp[,13] <- sapply(meta.samp[,13] ,as.numeric)

pdf(paste("Read_Length_v_Quality_marginals", suffix))
plot(ScatterHist(meta.samp, "Q_Score", "Read_Length",
            annot_size=10,
            title="Read Length v. Quality Score"))
dev.off()     
       
meta[,14] <- sapply(meta[,14] ,as.character)
meta[,14] <- sapply(meta[,14] ,as.numeric)
meta[,13] <- sapply(meta[,13] ,as.character)
meta[,13] <- sapply(meta[,13] ,as.numeric)
meta[,5] <- sapply(meta[,5] ,as.character)
meta[,5] <- sapply(meta[,5] ,as.numeric)
meta.sort.time=meta[order(as.numeric(as.character(meta$start_time))), ]

meta.t.mod=meta
sp=as.numeric(strsplit(suffix2, "_")[[1]][3])
meta.t.mod$start_time[meta.t.mod$run_id=="5dbf043921978e278537bae7b93dd39e58f2e963"]=meta.t.mod$start_time[meta.t.mod$run_id=="5dbf043921978e278537bae7b93dd39e58f2e963"]+407.747+sp

meta.t.mod$start_time[meta.t.mod$run_id=="5eb70612726a8e9c120c71d557ab18757da2fdb8"]=meta.t.mod$start_time[meta.t.mod$run_id=="5eb70612726a8e9c120c71d557ab18757da2fdb8"]+407.747+sp+24960.9+33480

meta.t.mod$start_time[meta.t.mod$run_id=="d107a32592fb48820492187597eec7acd795c319"]=meta.t.mod$start_time[meta.t.mod$run_id=="d107a32592fb48820492187597eec7acd795c319"]+407.747+sp+24960.9+33480+403.0805+sp

meta.t.mod$start_time[meta.t.mod$run_id=="f2171b80d1e729a4d137fccb31087ac33c1c913e"]=meta.t.mod$start_time[meta.t.mod$run_id=="f2171b80d1e729a4d137fccb31087ac33c1c913e"]+407.747+sp+24960.9+33480+403.0805+sp+10719.92+sp

meta.t.mod$start_time[meta.t.mod$run_id=="d19bf639e1fdf76d71d4ad3979979d2e7fe1ea83"]=meta.t.mod$start_time[meta.t.mod$run_id=="d19bf639e1fdf76d71d4ad3979979d2e7fe1ea83"]+407.747+sp+24960.9+33480+403.0805+sp+10719.92+sp+396.5625+sp

meta.t.mod=meta.t.mod[order(meta.t.mod$start_time), ]

step=dim(meta.t.mod)[1]/10000
step=1

i=1
meta.s=0
meta.s.m=matrix(nrow=round(dim(meta.t.mod)[1]/step)+1,ncol=3)
i2=i+step
a=1
while (i<dim(meta.t.mod)[1]){
    meta.s=meta.s+sum(meta.t.mod[i:i2, 13])
    tstamp= meta.t.mod[i, 5]
    tstamp=tstamp/3600
    numb=a*step
    meta.s.m[a,1]=tstamp
    meta.s.m[a,2]=as.numeric(meta.s)
    meta.s.m[a,3]=numb
    i2=i+step+1
    i=i+step+1
    a=a+1
}

sub=lm(meta.s.m[,2]~meta.s.m[,1])

pdf(file=paste("Cumulative_Read_Length_Time",suffix2, suffix))

plot(meta.s.m, xlab="Time (hours)", ylab="Cumulative Read Length", main="Cumulative Read Length over Time", sub=paste("Coefficient: ",format(signif(sub$coefficients[2],3),scientific=TRUE), "bases per hour, ", "Intercept: ",format(signif(sub$coefficients[1],3), scientific=TRUE), ", R sq.: ",signif(as.numeric(summary(sub)$r.squared), 3)))
 
abline(lm(meta.s.m[,2]~meta.s.m[,1]), col="red")

dev.off()

sub=lm(meta.s.m[,3]~meta.s.m[,1])

pdf(file=paste("Cumulative_Read_Count_Time", suffix2, suffix))

plot(meta.s.m[,1], meta.s.m[,3], xlab="Time (hours)", ylab="Cumulative Read Count", main="Cumulative Read Count over Time", sub=paste("Coefficient: ",format(signif(sub$coefficients[2],3),scientific=TRUE), "reads per hour, ", "Intercept: ",format(signif(sub$coefficients[1],3), scientific=TRUE), ", R sq.: ",signif(as.numeric(summary(sub)$r.squared), 3)))

abline(lm(meta.s.m[,3]~meta.s.m[,1]), col="red")

dev.off()