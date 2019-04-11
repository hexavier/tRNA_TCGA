library(ggplot2)
library(ggpubr)

##### TOP #####
## Load codon centric data
cu_TAIs = read.csv("results/TOPweighted_tAI_CUcentric.csv",row.names = 1)

dataset_CU = c()
types = sapply(colnames(cu_TAIs), function(x) strsplit(x,"\\.")[[1]][1])

for (sample in colnames(cu_TAIs)){
  dataset_temp = data.frame(row.names=seq(1,ncol(cu_TAIs)))
  dataset_temp$tai = as.numeric(cu_TAIs[,sample])
  dataset_temp$match = (types %in% strsplit(sample,"\\.")[[1]][1])
  dataset_temp$sample_exp = colnames(cu_TAIs)
  dataset_temp$catype_exp = types
  dataset_CU = rbind(dataset_CU, dataset_temp)
}

## Load trna centric data
trna_TAIs = read.csv("results/TOPweighted_tAI_trnacentric.csv",row.names = 1)

dataset_trna = c()
types = sapply(colnames(trna_TAIs), function(x) strsplit(x,"\\.")[[1]][1])

for (sample in colnames(trna_TAIs)){
  dataset_temp = data.frame(row.names=seq(1,ncol(trna_TAIs)))
  dataset_temp$tai = as.numeric(trna_TAIs[,sample])
  dataset_temp$match = (types %in% strsplit(sample,"\\.")[[1]][1])
  dataset_temp$sample_tai = sample
  dataset_temp$catype_tai = strsplit(sample,"\\.")[[1]][1]
  dataset_trna = rbind(dataset_trna, dataset_temp)
}

## Build dataset for paired test
dataset_paired = c()
for (sample in colnames(cu_TAIs)){
  match = (types %in% strsplit(sample,"\\.")[[1]][1])
  dataset_temp = data.frame(row.names=seq(4))
  # Record 4 values: matched and unmatched trnas and cu.
  dataset_temp$tai = c(median(as.numeric(cu_TAIs[sample,!match]),na.rm = T),median(as.numeric(cu_TAIs[sample,match]),na.rm = T),
                       median(as.numeric(trna_TAIs[!match,sample]),na.rm = T),median(as.numeric(trna_TAIs[match,sample]),na.rm = T))
  dataset_temp$match = c(FALSE,TRUE,FALSE,TRUE)
  dataset_temp$comparison = c("exp","exp","tai","tai")
  dataset_temp$sample = sample
  dataset_temp$catype = strsplit(sample,"\\.")[[1]][1]
  dataset_paired = rbind(dataset_paired, dataset_temp)
}
compare_means(tai ~ match,data=dataset_paired,group.by = "comparison",method="wilcox.test",alternative = "greater",paired=T)

pdf("plots/tAI_TOP_paired.pdf",height=8,width = 12)
print(ggpaired(dataset_paired, x = "match", y = "tai",color = "match", id="sample",facet.by = "comparison", line.color="gray", xlab="Match", ylab="tAI", title="tAI between tissues") + stat_compare_means(label = "p.format", method.args = list(alternative = "greater"), paired = TRUE))
print(ggpaired(dataset_paired[dataset_paired$comparison=="exp",], x = "match", y = "tai",color = "match", id="sample",facet.by = "catype", line.color="gray", xlab="Match", ylab="tAI", title="tAI between tissues") + stat_compare_means(label = "p.format", method.args = list(alternative = "greater"), paired = TRUE))
print(ggpaired(dataset_paired[dataset_paired$comparison=="tai",], x = "match", y = "tai",color = "match", id="sample",facet.by = "catype", line.color="gray", xlab="Match", ylab="tAI", title="tAI between tissues") + stat_compare_means(label = "p.format", method.args = list(alternative = "greater"), paired = TRUE))
dev.off()

pdf("plots/tAI_TOP_unpaired.pdf",height=5,width = 15)
print(ggplot(dataset_trna, aes(x=catype_tai, y=tai, fill=match)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="tAI between tissues tRNA",x="Tissue", y = "tAI") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())
print(ggplot(dataset_CU, aes(x=catype_exp, y=tai, fill=match)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="tAI between tissues CU",x="Tissue", y = "tAI") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())
dev.off()

##### BOTTOM #####
## Load codon centric data
cu_TAIs = read.csv("results/BOTTOMweighted_tAI_CUcentric.csv",row.names = 1)

dataset_CU = c()
types = sapply(colnames(cu_TAIs), function(x) strsplit(x,"\\.")[[1]][1])

for (sample in colnames(cu_TAIs)){
  dataset_temp = data.frame(row.names=seq(1,ncol(cu_TAIs)))
  dataset_temp$tai = as.numeric(cu_TAIs[,sample])
  dataset_temp$match = (types %in% strsplit(sample,"\\.")[[1]][1])
  dataset_temp$sample_exp = colnames(cu_TAIs)
  dataset_temp$catype_exp = types
  dataset_CU = rbind(dataset_CU, dataset_temp)
}

## Load trna centric data
trna_TAIs = read.csv("results/BOTTOMweighted_tAI_trnacentric.csv",row.names = 1)

dataset_trna = c()
types = sapply(colnames(trna_TAIs), function(x) strsplit(x,"\\.")[[1]][1])

for (sample in colnames(trna_TAIs)){
  dataset_temp = data.frame(row.names=seq(1,ncol(trna_TAIs)))
  dataset_temp$tai = as.numeric(trna_TAIs[,sample])
  dataset_temp$match = (types %in% strsplit(sample,"\\.")[[1]][1])
  dataset_temp$sample_tai = sample
  dataset_temp$catype_tai = strsplit(sample,"\\.")[[1]][1]
  dataset_trna = rbind(dataset_trna, dataset_temp)
}

## Build dataset for paired test
dataset_paired = c()
for (sample in colnames(cu_TAIs)){
  match = (types %in% strsplit(sample,"\\.")[[1]][1])
  dataset_temp = data.frame(row.names=seq(4))
  # Record 4 values: matched and unmatched trnas and cu.
  dataset_temp$tai = c(median(as.numeric(cu_TAIs[sample,!match]),na.rm = T),median(as.numeric(cu_TAIs[sample,match]),na.rm = T),
                       median(as.numeric(trna_TAIs[!match,sample]),na.rm = T),median(as.numeric(trna_TAIs[match,sample]),na.rm = T))
  dataset_temp$match = c(FALSE,TRUE,FALSE,TRUE)
  dataset_temp$comparison = c("exp","exp","tai","tai")
  dataset_temp$sample = sample
  dataset_temp$catype = strsplit(sample,"\\.")[[1]][1]
  dataset_paired = rbind(dataset_paired, dataset_temp)
}
compare_means(tai ~ match,data=dataset_paired,group.by = "comparison",method="wilcox.test",alternative = "greater",paired=T)

pdf("plots/tAI_BOTTOM_paired.pdf",height=8,width = 12)
print(ggpaired(dataset_paired, x = "match", y = "tai",color = "match", id="sample",facet.by = "comparison", line.color="gray", xlab="Match", ylab="tAI", title="tAI between tissues") + stat_compare_means(label = "p.format", method.args = list(alternative = "greater"), paired = TRUE))
print(ggpaired(dataset_paired[dataset_paired$comparison=="exp",], x = "match", y = "tai",color = "match", id="sample",facet.by = "catype", line.color="gray", xlab="Match", ylab="tAI", title="tAI between tissues") + stat_compare_means(label = "p.format", method.args = list(alternative = "greater"), paired = TRUE))
print(ggpaired(dataset_paired[dataset_paired$comparison=="tai",], x = "match", y = "tai",color = "match", id="sample",facet.by = "catype", line.color="gray", xlab="Match", ylab="tAI", title="tAI between tissues") + stat_compare_means(label = "p.format", method.args = list(alternative = "greater"), paired = TRUE))
dev.off()

pdf("plots/tAI_BOTTOM_unpaired.pdf",height=5,width = 15)
print(ggplot(dataset_trna, aes(x=catype_tai, y=tai, fill=match)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="tAI between tissues tRNA",x="Tissue", y = "tAI") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())
print(ggplot(dataset_CU, aes(x=catype_exp, y=tai, fill=match)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="tAI between tissues CU",x="Tissue", y = "tAI") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())
dev.off()
