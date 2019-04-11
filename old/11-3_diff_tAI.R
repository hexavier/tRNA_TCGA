library(ggplot2)
library(ggpubr)

## Load data
TAIs = read.csv("results/weighted_tAI.csv",row.names = 1)
modTAIs = read.csv("results/weighted_modif_tAI.csv",row.names = 1)

##### PLOT VERSUS TYPES #####
dataset = c()
types = sapply(colnames(TAIs), function(x) strsplit(x,"\\.")[[1]][1])

for (sample in colnames(TAIs)){
  dataset_temp = data.frame(row.names=seq(1,2))
  dataset_temp$tai = as.numeric(c(TAIs[sample,sample],modTAIs[sample,sample]))
  dataset_temp$modif = c(FALSE,TRUE)
  dataset_temp$sample = sample
  dataset_temp$catype = strsplit(sample,"\\.")[[1]][1]
  dataset = rbind(dataset, dataset_temp)
}

compare_means(tai ~ modif,data=dataset,group.by = "catype",method="wilcox.test",alternative = "greater",paired=T)

pdf("plots/tAI_modVSnomod_paired.pdf",height=8,width = 12)
print(ggpaired(dataset, x = "modif", y = "tai",color = "modif", id="sample", line.color="gray", xlab="Match", ylab="tAI", title="tAI between tissues") + stat_compare_means(label = "p.format", method.args = list(alternative = "greater"), paired = TRUE))
print(ggpaired(dataset, x = "modif", y = "tai",color = "modif", id="sample",facet.by = "catype", line.color="gray", xlab="Match", ylab="tAI", title="tAI between tissues") + stat_compare_means(label = "p.format", method.args = list(alternative = "greater"), paired = TRUE))
dev.off()

## Analyze tAI matched vs mismatch in modifies dataset
dataset = c()
types = sapply(colnames(modTAIs), function(x) strsplit(x,"\\.")[[1]][1])

for (sample in colnames(modTAIs)){
  dataset_temp = data.frame(row.names=seq(1,ncol(modTAIs)))
  dataset_temp$tai = as.numeric(modTAIs[sample,])
  dataset_temp$match = (types %in% strsplit(sample,"\\.")[[1]][1])
  dataset_temp$sample_exp = sample
  dataset_temp$sample_tai = colnames(modTAIs)
  dataset_temp$catype_exp = strsplit(sample,"\\.")[[1]][1]
  dataset_temp$catype_tai = types
  dataset = rbind(dataset, dataset_temp)
}


pdf("plots/tAI_modif_unpaired.pdf",height=5,width = 15)
print(ggplot(dataset, aes(x=catype_tai, y=tai, fill=match)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="tAI between tissues tRNAs",x="Tissue", y = "tAI") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())
print(ggplot(dataset, aes(x=catype_exp, y=tai, fill=match)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="tAI between tissues CU",x="Tissue", y = "tAI") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())
dev.off()
