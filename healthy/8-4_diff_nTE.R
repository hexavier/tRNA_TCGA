library(ggplot2)
library(ggpubr)

# Load nTE data
load(file="/home/xhernandez/Documents/tAI_genomic/AAcTE_CUgenomic_sqrt.rda")
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# Extract matched samples
sums = apply(nTE,2,rowSums)
# Rename samples based on clustering
colnames(sums) = gsub("KICH","kidney",colnames(sums))
colnames(sums) = gsub("KIRP","kidney",colnames(sums))
colnames(sums) = gsub("KIRC","kidney",colnames(sums))
colnames(sums) = gsub("LUAD","lung",colnames(sums))
colnames(sums) = gsub("LUSC","lung",colnames(sums))
colnames(sums) = gsub("COAD","colorectal",colnames(sums))
colnames(sums) = gsub("READ","colorectal",colnames(sums))
colnames(sums) = gsub("CHOL","liver",colnames(sums))
colnames(sums) = gsub("LIHC","liver",colnames(sums))
colnames(sums) = gsub("UCEC","uterus",colnames(sums))
colnames(sums) = gsub("CESC","uterus",colnames(sums))

rownames(sums) = gsub("KICH","kidney",rownames(sums))
rownames(sums) = gsub("KIRP","kidney",rownames(sums))
rownames(sums) = gsub("KIRC","kidney",rownames(sums))
rownames(sums) = gsub("LUAD","lung",rownames(sums))
rownames(sums) = gsub("LUSC","lung",rownames(sums))
rownames(sums) = gsub("COAD","colorectal",rownames(sums))
rownames(sums) = gsub("READ","colorectal",rownames(sums))
rownames(sums) = gsub("CHOL","liver",rownames(sums))
rownames(sums) = gsub("LIHC","liver",rownames(sums))
rownames(sums) = gsub("UCEC","uterus",rownames(sums))
rownames(sums) = gsub("CESC","uterus",rownames(sums))

##### PLOT VERSUS TYPES #####
dataset = c()
types = sapply(colnames(sums), function(x) strsplit(x,"\\.")[[1]][1])

for (sample in colnames(sums)){
  dataset_temp = data.frame(row.names=seq(1,ncol(sums)))
  dataset_temp$tai = as.numeric(sums[sample,])
  dataset_temp$match = (types %in% strsplit(sample,"\\.")[[1]][1])
  dataset_temp$sample_exp = sample
  dataset_temp$sample_tai = colnames(sums)
  dataset_temp$catype_exp = strsplit(sample,"\\.")[[1]][1]
  dataset_temp$catype_tai = types
  dataset = rbind(dataset, dataset_temp)
}

compare_means(tai ~ match,data=dataset,group.by = "catype_tai",method="wilcox.test",alternative = "less")
compare_means(tai ~ match,data=dataset,group.by = "catype_exp",method="wilcox.test",alternative = "less")

pdf("plots/AAcTE_codons_unpaired.pdf",height=5,width = 15)
print(ggplot(dataset, aes(x=catype_tai, y=tai, fill=match)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="nTE between tissues tRNAs",x="Tissue", y = "Sum TE") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())
print(ggplot(dataset, aes(x=catype_exp, y=tai, fill=match)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="nTE between tissues CU",x="Tissue", y = "Sum TE") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())

##### PLOT MISMATCH TOGETHER #####
compare_means(tai ~ match,data=dataset,method="wilcox.test",alternative = "less")

print(ggplot(dataset, aes(x=match, y=tai)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="nTE between tissues CU",x="Match", y = "Sum TE") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())
dev.off()
