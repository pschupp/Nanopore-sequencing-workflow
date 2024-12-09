genes =read.csv("Oldham_fidelity_4_celltypes_genes.csv")
genes=genes[,-1]
expr=read.table("nano_24_sec_seq_expression_matrix_log2_norm_k5.tsv", sep="\t", header=TRUE)



for (types in seq(1,4)){
print(colnames(genes)[types])
print(length(which(genes[,types] %in% rownames(expr))))}

[1] "Astrocyte"
[1] 32
[1] "Microglia"
[1] 11
[1] "Neuron"
[1] 25
[1] "Oligodendrocyte"
[1] 30

for (types in seq(1,4)){
which(rownames(expr) %in% genes[,types])
expr.working=apply(expr[which(rownames(expr) %in% genes[,types]), ], 1, sum)
print(colnames(genes)[types])
print(expr.working)
}
[1] "Astrocyte"
  ATP1A2   SLC1A3    TIMP3     MLC1     PON2   SLC1A2    ALDH2  ALDH6A1 
     440      133      930       26     1100     1259       13     1207 
    SOX9  ATP13A4     APOE   NOTCH2     LRP4      AGT    EDNRB   BMPR1B 
      52       44       65        3       91        6       48     4529 
 TP53BP2    MYO10    NTRK2     GJA1   PDLIM5   ETNPPL    RAB31   AHCYL1 
       7       37      358       17       35       15       59      169 
   S1PR1     AQP4 CDC42EP4     SOX2 SLC25A18  METTL7A    SPON1 ARHGEF26 
    1901        7       17       81       11      185        7       27 
[1] "Microglia"
RNASET2      C3  RHBDF2  HAVCR2   ITGB2     SYK   PTAFR   GPR34   CSF1R   ADAP2 
     12       5       6       8       9      11       9       8      16      10 
HLA-DRA 
     16 
[1] "Neuron"
   CELF2  GUCY1B3     ADD2    MAPK1   TRIM37     CAP2      NNT   PPP2CA 
      36        9      494       29       17       17        7       25 
   RNF11     NAPB    RHOT1   NDFIP1   SNAP25    MYH10    G3BP2     ANK2 
       5       14       15       18       22       25       33      167 
ATP6V1B2    CADPS   UBE2V2    NRXN1    SCN8A   FAM49A     OPA1   BHLHB9 
      21      101       70       54        8        9      358       51 
     NSF 
       8 
[1] "Oligodendrocyte"
   ERBB3  FAM107B  SLC44A1    PSEN1   CLDND1       TF     RFFL     FA2H 
       8       15      894       42       16       19       14        4 
   PMP22     CA14     PLP1    HSPA2   RNASE1   SLAIN1     NPC1  PIP4K2A 
       6        4       43       13        7        8       84     6122 
C10orf90  TMEM144    USP54      PXK   SH3TC2    GPR37   CARNS1      CNP 
     176      507       56        8       29        8      375     1914 
    UGT8  PLA2G16     GLDN   ZDHHC9  TMEM63A     NDE1 
       7       70      393        4       40       13 
       
       
for (types in seq(1,4)){
which(rownames(expr) %in% genes[,types])
expr.working=apply(expr[which(rownames(expr) %in% genes[,types]), ], 1, sum)
print(colnames(genes)[types])
print(sum(expr.working))
}

[1] "Astrocyte"
[1] 12879
[1] "Microglia"
[1] 110
[1] "Neuron"
[1] 1613
[1] "Oligodendrocyte"
[1] 10899