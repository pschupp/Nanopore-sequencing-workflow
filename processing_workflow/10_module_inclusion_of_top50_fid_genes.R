fid=read.csv("~/Documents/oldham_gene_interest/top50_fid_genes.csv")
kme=read.csv("kME_table_03-16-44.csv")

modules=sort(unique(kme$TopModPosFDR_0.000635)) #change fdr pval for all mentions if using different file
modules=sapply(modules, as.character)
out.dat=as.data.frame(matrix(nrow=length(modules), ncol=5))
out.dat[,1]=modules
colnames(out.dat)=c("Modules", "Astrocytes", "Microglia", "Neurons", "Oligodendrocytes")

m=1
for (module in out.dat$Modules){        #this snippet of code get all the top 50 fidelty genes that are present in each of the 4 gene catagories
    kme.w=kme[which(kme$TopModPosFDR_0.000635==module),c(1, 2,3,4,5,6)]
    c=2
    for (type in colnames(fid)){
        num=as.character(kme.w[(which(kme.w[,1] %in% fid[,c-1])),1])
        out.dat[m,c]=paste(num, collapse='_')
        c=c+1
    }
    m=m+1
} 

exc=vector('list', 4) #this snippet of code gets all the genes that are not in any module for each of the 4 cell types
c=1
list=kme[which(mod != "NA"),1]

for (type in colnames(fid)){
    num=as.character(fid[(which(!(fid[,c] %in% list))),c])
    exc[[as.character(type)]]=num
    c=c+1
}

#want to compare avg expression values for top200 fid genes that are in modules and the remaining ones which are not
mod=as.character(kme[,6])
genes.mod=kme[which(mod != "NA"),1]

out=c()
for (type in colnames(fid)){
    inc.genes=fid[[type]][which(fid[[type]] %in% genes.mod)]
    out[[as.character(type)]]=inc.genes
}

expr=read.table("~/Documents/24_sec_seq/root/09_expression_matrix/GRCh38/RUVg_count_data_5_top_ERCC/nano_24_sec_seq_expression_matrix_log2_norm_k5.tsv", sep="\t", header=TRUE)

expr.comp=c()
x=1
for (types in out){
    expr.w=expr[which(rownames(expr) %in% types),]
    expr.w=apply(expr.w, 1, sum)
    expr.comp[[colnames(fid)[x]]]=unlist(expr.w)
    x=x+1
}
x=1
for (types in out){
    expr.w=expr[which(rownames(expr) %in% exc[[colnames(fid)[x]]]) , ]
    expr.w=apply(expr.w, 1, sum)
    expr.w=expr.w[sample(seq(1:length(expr.w)),length(expr.comp[[x]]))]
    expr.comp[[paste(colnames(fid)[x], "un_inc")]]=unlist(expr.w)
    x=x+1
}

pdf("expression_fid_mod_adj.pdf",width=15)
boxplot(expr.comp, xlab="Top 50 Fidelity Genes for Each Cell Type Included in Modules / Top 50 Fidelity Genes for Each Cell Type Not Included in Modules", ylab="Counts (raw)", main="Expression for the top 50 fidelity genes for each cell type included or not included in modules, not included gene number downsampled to match included")
for (x in seq(1,8)){
    points(rep(x,length(expr.comp[[x]])), expr.comp[[x]], col='red')}
dev.off()


expr.comp=c()
x=1
for (types in out){
    expr.w=expr[which(rownames(expr) %in% types),]
    expr.w=apply(expr.w, 1, sum)
    expr.comp[[colnames(fid)[x]]]=unlist(expr.w)
    x=x+1
}
x=1
for (types in out){
    expr.w=expr[which(rownames(expr) %in% exc[[colnames(fid)[x]]]) , ]
    expr.w=apply(expr.w, 1, sum)
   # expr.w=expr.w[sample(seq(1:length(expr.w)),length(expr.comp[[x]]))]
    expr.comp[[paste(colnames(fid)[x], "un_inc")]]=unlist(expr.w)
    x=x+1
}

pdf("expression_fid_mod_unadj.pdf",width=15)
boxplot(expr.comp, xlab="Top 50 Fidelity Genes for Each Cell Type Included in Modules / Top 50 Fidelity Genes for Each Cell Type Not Included in Modules", ylab="Counts (raw)",  main="Expression for the top 50 fidelity genes for each cell type included or not included in modules")
for (x in seq(1,8)){
    points(rep(x,length(expr.comp[[x]])), expr.comp[[x]], col='red')}
dev.off()