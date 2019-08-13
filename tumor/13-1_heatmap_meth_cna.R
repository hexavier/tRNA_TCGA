library(gplots)
library(ggplot2)

#%% Load data
cancer_types = c(BRCA="BRCA",PRAD="PRAD",KICH="KICH",KIRC="KIRC",KIRP="KIRP",LUAD="LUAD",LUSC="LUSC",HNSC="HNSC",UCEC="UCEC",CESC="CESC",
                 LIHC="LIHC", CHOL="CHOL",THCA="THCA",COAD="COAD",READ="READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG")
#values = read.csv("results/differentialMeth_TvsH_twosided.csv",row.names = 1)
#values = read.csv("results/differentialCNA_TvsH_twosided.csv",row.names = 1)
#values = read.csv("results/differentialMeanMeth_TvsH_twosided.csv",row.names = 1)
values = read.csv("results/differentialMeanCNA_TvsH_twosided.csv",row.names = 1)

#aadata = values[substr(values$isoacceptor,1,3)=="Gly",]
aadata = values[substr(values$gene,1,3)=="Gly",]
genes = unique(aadata$gene)

# Initialize structures
diffexp_mean = matrix(nrow = length(genes), ncol = length(cancer_types))
rownames(diffexp_mean) = genes; colnames(diffexp_mean) = names(cancer_types)
is_sig = matrix(nrow = length(genes),ncol = length(cancer_types))
rownames(is_sig) = genes; colnames(is_sig) = names(cancer_types)

for (type in names(cancer_types)){
  # Add data
  diffexp_mean[,type] = sapply(genes, function(x) aadata[(aadata$catype==type)&(aadata$gene==x),"FC"])
  is_sig[,type] = sapply(genes, function(x) aadata[(aadata$catype==type)&(aadata$gene==x),"pval_adj"])
}

# Keep only significant exons
is_sig = is_sig[rowSums(!is.na(diffexp_mean))>0,]
diffexp_mean = diffexp_mean[rowSums(!is.na(diffexp_mean))>0,]
diffexp_mean[(is_sig>=0.05)|(is.na(is_sig))] = NA

# Heatmap
breaks = c(min(diffexp_mean,na.rm=T),seq(-1,1,length.out=254),max(diffexp_mean,na.rm=T))
diffexp_mean = diffexp_mean[order((rowSums(diffexp_mean<0,na.rm=T)-rowSums(diffexp_mean>0,na.rm=T)),decreasing = T),order(colSums(is.na(diffexp_mean)))]
heatmap.2(diffexp_mean,
          Rowv=F, Colv=F, col = bluered(255), symm=T, margins=c(12,15), trace="none",
          na.color="grey",breaks = breaks)
