library(ggplot2)
library(ggpubr)

# Load nTE data
load(file="/home/xhernandez/Documents/tAI_genomic/AAcTE_CUgenomic_sqrt.rda")

# Rename samples based on clustering
colnames(nTE) = gsub("KICH","C1",colnames(nTE))
colnames(nTE) = gsub("KIRP","C1",colnames(nTE))
colnames(nTE) = gsub("KIRC","C1",colnames(nTE))
colnames(nTE) = gsub("LUAD","C2",colnames(nTE))
colnames(nTE) = gsub("LUSC","C2",colnames(nTE))
colnames(nTE) = gsub("COAD","C3",colnames(nTE))
colnames(nTE) = gsub("READ","C3",colnames(nTE))
colnames(nTE) = gsub("CHOL","C1",colnames(nTE))
colnames(nTE) = gsub("LIHC","C1",colnames(nTE))
colnames(nTE) = gsub("UCEC","C3",colnames(nTE))
colnames(nTE) = gsub("CESC","C3",colnames(nTE))
colnames(nTE) = gsub("GBM","C1",colnames(nTE))
colnames(nTE) = gsub("PCPG","C1",colnames(nTE))
colnames(nTE) = gsub("THCA","C2",colnames(nTE))
colnames(nTE) = gsub("PRAD","C3",colnames(nTE))
colnames(nTE) = gsub("BRCA","C2",colnames(nTE))
colnames(nTE) = gsub("BLCA","C3",colnames(nTE))
colnames(nTE) = gsub("ESCA","C1",colnames(nTE))
colnames(nTE) = gsub("PAAD","C3",colnames(nTE))
colnames(nTE) = gsub("STAD","C3",colnames(nTE))
colnames(nTE) = gsub("HNSC","C2",colnames(nTE))
colnames(nTE) = gsub("THYM","C2",colnames(nTE))
colnames(nTE) = gsub("SKCM","C3",colnames(nTE))

rownames(nTE) = gsub("KICH","C1",rownames(nTE))
rownames(nTE) = gsub("KIRP","C1",rownames(nTE))
rownames(nTE) = gsub("KIRC","C1",rownames(nTE))
rownames(nTE) = gsub("LUAD","C2",rownames(nTE))
rownames(nTE) = gsub("LUSC","C2",rownames(nTE))
rownames(nTE) = gsub("COAD","C3",rownames(nTE))
rownames(nTE) = gsub("READ","C3",rownames(nTE))
rownames(nTE) = gsub("CHOL","C1",rownames(nTE))
rownames(nTE) = gsub("LIHC","C1",rownames(nTE))
rownames(nTE) = gsub("UCEC","C3",rownames(nTE))
rownames(nTE) = gsub("CESC","C3",rownames(nTE))
rownames(nTE) = gsub("GBM","C1",rownames(nTE))
rownames(nTE) = gsub("PCPG","C1",rownames(nTE))
rownames(nTE) = gsub("THCA","C2",rownames(nTE))
rownames(nTE) = gsub("PRAD","C3",rownames(nTE))
rownames(nTE) = gsub("BRCA","C2",rownames(nTE))
rownames(nTE) = gsub("BLCA","C3",rownames(nTE))
rownames(nTE) = gsub("ESCA","C1",rownames(nTE))
rownames(nTE) = gsub("PAAD","C3",rownames(nTE))
rownames(nTE) = gsub("STAD","C3",rownames(nTE))
rownames(nTE) = gsub("HNSC","C2",rownames(nTE))
rownames(nTE) = gsub("THYM","C2",rownames(nTE))
rownames(nTE) = gsub("SKCM","C3",rownames(nTE))

##### PLOT VERSUS TYPES #####
sums = apply(nTE,2,rowSums)
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

pdf("plots/AAcTE_bycluster.pdf",height=5,width = 5)
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
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "less")) +
        theme_classic())
dev.off()

##### PLOT VERSUS TYPES #####
# Detect types and matching tissues
types = sapply(colnames(nTE), function(x) strsplit(x,"\\.")[[1]][1])
matched = as.logical(sapply(types,function(s)(types %in% s)))

# Build dataset
dataset = data.frame(row.names = 1:((ncol(nTE)^2)*60))
dataset$tai = as.numeric(nTE)
dataset$match = rep(matched,60)
dataset$codon = as.character(sapply(dimnames(nTE)[[3]], rep, (ncol(nTE)^2)))

tests = compare_means(tai ~ match,data=dataset,group.by = "codon",method="wilcox.test",alternative = "less")

# Plots
pdf("plots/AAcTE_bycluster_bycod.pdf",height=5,width = 5)
for (c in dimnames(nTE)[[3]]){
  # Retrieve data
  dataset = data.frame(row.names = 1:(ncol(nTE)^2))
  dataset$tai = as.numeric(nTE[,,c])
  dataset$match = matched
  dataset$codon = c
  # Plot
  print(ggplot(dataset, aes(x=codon, y=tai, fill=match)) +  
          geom_violin(position=position_dodge(1)) +
          labs(title="TE between codons",x="Codon", y = "TE") + 
          stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
          stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "less")) +
          theme_classic())
}
dev.off()
