library(ggplot2)
library(ggpubr)

pdf("plots/TS_AAcTE_codons_unpaired.pdf",height=5,width = 15)

##### SUPPLY #####
# Load nTE data
load(file="/home/xhernandez/Documents/tAI_genomic/AAcTE_TSCUgenomic_sqrt_SUPcentric.rda")
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# Extract matched samples
sums = apply(nTE,2,rowSums)

##### PLOT VERSUS TYPES #####
dataset = c()
types = sapply(colnames(sums), function(x) strsplit(x,"\\.")[[1]][1])

for (sample in colnames(sums)){
  dataset_temp = data.frame(row.names=seq(1,ncol(sums)))
  dataset_temp$tai = as.numeric(sums[sample,])
  dataset_temp$match = (types %in% strsplit(sample,"\\.")[[1]][1])
  dataset_temp$sample_tai = colnames(sums)
  dataset_temp$catype_tai = types
  dataset = rbind(dataset, dataset_temp)
}

compare_means(tai ~ match,data=dataset,group.by = "catype_tai",method="wilcox.test",alternative = "less")

print(ggplot(dataset, aes(x=catype_tai, y=tai, fill=match)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="nTE between tissues SUPPLY",x="Tissue", y = "Sum TE") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())

compare_means(tai ~ match,data=dataset,method="wilcox.test",alternative = "less")

print(ggplot(dataset, aes(x=match, y=tai)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="nTE SUPPLY",x="Match", y = "Sum TE") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())

##### DEMAND #####
# Load nTE data
load(file="/home/xhernandez/Documents/tAI_genomic/AAcTE_TSCUgenomic_sqrt_DEMcentric.rda")
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# Extract matched samples
sums = apply(nTE,2,rowSums)

##### PLOT VERSUS TYPES #####
dataset = c()
types = sapply(colnames(sums), function(x) strsplit(x,"\\.")[[1]][1])

for (sample in colnames(sums)){
  dataset_temp = data.frame(row.names=seq(1,ncol(sums)))
  dataset_temp$tai = as.numeric(sums[sample,])
  dataset_temp$match = (types %in% strsplit(sample,"\\.")[[1]][1])
  dataset_temp$sample_exp = sample
  dataset_temp$catype_exp = strsplit(sample,"\\.")[[1]][1]
  dataset = rbind(dataset, dataset_temp)
}

compare_means(tai ~ match,data=dataset,group.by = "catype_exp",method="wilcox.test",alternative = "less")

print(ggplot(dataset, aes(x=catype_exp, y=tai, fill=match)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="nTE between tissues DEMAND",x="Tissue", y = "Sum TE") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())

compare_means(tai ~ match,data=dataset,method="wilcox.test",alternative = "less")

print(ggplot(dataset, aes(x=match, y=tai)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="nTE DEMAND",x="Tissue", y = "Sum TE") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())

dev.off()

