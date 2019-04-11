library(ggplot2)
library(ggpubr)

##### ANALYZE TOP EXPRESSED GENES #####
## Load data
TAIs = read.csv("results/CUbybins/tAI_relcodons/tAI_relcodons_bin0.csv",row.names = 1)

# PLOT VERSUS TYPES
dataset = c()
types = sapply(colnames(TAIs), function(x) strsplit(x,"\\.")[[1]][1])

for (sample in colnames(TAIs)){
  dataset_temp = data.frame(row.names=seq(1,ncol(TAIs)))
  dataset_temp$tai = as.numeric(TAIs[sample,])
  dataset_temp$match = (types %in% strsplit(sample,"\\.")[[1]][1])
  dataset_temp$sample_trnas = sample
  dataset_temp$sample_codon = colnames(TAIs)
  dataset_temp$catype_trnas = strsplit(sample,"\\.")[[1]][1]
  dataset_temp$catype_codon = types
  dataset = rbind(dataset, dataset_temp)
}

pdf("plots/tAI_relcodons_TOP_bytype.pdf",height=5,width = 15)
print(ggplot(dataset, aes(x=catype_trnas, y=tai, fill=match)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="tAI between tissues tRNAs",x="Tissue", y = "tAI") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())
print(ggplot(dataset, aes(x=catype_codon, y=tai, fill=match)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="tAI between tissues CU",x="Tissue", y = "tAI") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",method.args = list(alternative = "greater")) +
        theme_classic())
dev.off()

# PLOT INDIVIDUALLY
dataset = c()
types = sapply(colnames(TAIs), function(x) strsplit(x,"\\.")[[1]][1])

for (sample in colnames(TAIs)){
  dataset_temp = data.frame(row.names=sample)
  match = as.numeric(TAIs[sample,sample])
  dataset_temp$taimatch = match
  nomatch = as.numeric(TAIs[sample,!(colnames(TAIs) %in% sample)])
  dataset_temp$mediantainomatch = median(nomatch,na.rm = T)
  dataset_temp$randtest = sum(nomatch >= match,na.rm=T)/length(nomatch)
  dataset_temp$catype = strsplit(sample,"\\.")[[1]][1]
  dataset = rbind(dataset, dataset_temp)
}

write.csv(dataset,"results/diff_tAI_rel_TOP.csv")

##### ANALYZE tAI ACROSS BINS #####

# Init structure
tai_bins = data.frame(row.names = colnames(TAIs))
# Load data
bins = 20
for (n in 1:bins){
  taibin = read.csv(sprintf("results/CUbybins/tAI_relcodons/tAI_relcodons_bin%s.csv",(n-1)),row.names = 1)
  tai_bins[,n] = diag(as.matrix(taibin))
}
corcoef = apply(tai_bins[,1:12],1,function(x) cor(1:12,as.numeric(x)))
