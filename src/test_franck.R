rm(list=ls())
library(pscl)
library(ZIM)
library(glmmADMB)
library(parallel)
library(scRNAtools)
devtools::load_all("../scRNAtools/", reset = T)



#############################################################
X  = read.table("data/f_test/merged_peaks_counts.tsv",h=T)
peakid = X[,c(1:2)]
X  = X[,-c(1:2)]

# plot(density(log(apply(X,1,mean)))); rug(log(apply(X,1,mean)))
v = log(apply(X,1,mean))
X = X[(v<=exp(8)),]
peakid = peakid[(v<=exp(8)),]


hh = colnames(X)
o  = order(hh)
X  = X[,o]
hh = hh[o]
clonelist = c("P9855_01","P9855_02","P9855_03","P9855_04","P9855_05a","P9855_05b","P9855_07","P9855_08","P9855_09","P9855_10","P9855_11","P9855_13","P9855_14","P9855_15","P9855_16","P9855_17","P9855_18","P9855_19","P9855_20","P9855_21","P9855_22","P9855_23","P9855_24")
infos = data.frame(
    id = hh,
    clonality= rep(clonelist,c(2,1,1,1,1,1,1,2,2,1,1,1,1,2,1,1,1,1,1,1,2,2,1)),
    rep = c(1,2,1,1,1,1,1,1,1,2,1,2,1,1,1,1,1,2,1,1,1,1,1,1,1,2,1,2,1),
    clonecov = apply(X,2,mean)
    )
counts = matrix(t(X),nrow = dim(X)[2], ncol = dim(X)[1], dimnames = list(hh,peakid[,1]))
counts = counts[c(1,2,9,10,11,12,17,18,25,26,27,28),1:10]

library(xlsx)
A   = read.xlsx("data/f_test/P9855_ATAC_RNA_Worksheet_April4_2018n.xlsx",1)[,1:11]
colnames(A) = c("Sample.ID","Donor","Condition","HLA.Type","Batch","Clone.ID","Cell.Number","In_Vivo_Clone_Number","clonality","founder_pheno","protein_pheno")

tt        = apply(table(A$clonality,A$founder_pheno),1,which.max)
clonetype = data.frame(clonetype = sapply(tt,FUN=function(x){c("EFF","INT","MEM","none")[x]}))
clonetype$clonality  = row.names(clonetype)
row.names(clonetype) = NULL
infos = merge(infos,clonetype,by="clonality")

replicates = data.frame(clonality = names(which(table(infos$clonality)==2)))
infos = merge(infos,replicates,by="clonality")
infos <- infos[, c(2,1,3,4,5)]
scd = scdata$new(infos,counts)
dim(scd$getfeatures)
dim(scd$getcounts)
scd$getfeatures
scd$getcounts

system("mkdir -p results/f_test/")

system("rm -R results/f_test/")
out1 = DEA(scd,
    "y ~ 1",
    "y ~ clonetype",
    b_cells=rep(TRUE, scd$getncells),
    zi_threshold = 0.9,
    cpus = 1,
    v = T,
    folder_name = "results/f_test/"
  )
out1
hist(scd$getgene( "macs2_peak_5" ))


system("rm -R results/f_test/")
out2 = DEA(scd$select(genes = c("macs2_peak_5", "macs2_peak_6")),
    "y ~ 1 + offset(log(clonecov))",
    "y ~ clonality + offset(log(clonecov))",
    b_cells=rep(TRUE, scd$getncells),
    zi_threshold = 0.9,
    cpus = 1,
    v = T,
    folder_name = "results/f_test/",
    continuous = c("clonecov")
  )
out2

devtools::load_all("../scRNAtools/", reset = T)
system("rm -R results/f_test/")
out2 = DEA(scd,
    "y ~ 1 + offset(log(clonecov))",
    "y ~ clonality + offset(log(clonecov))",
    b_cells=rep(TRUE, scd$getncells),
    zi_threshold = 0.9,
    cpus = 1,
    v = T,
    folder_name = "results/f_test/",
    continuous = c("clonecov")
  )
traceback()

out2
out1
