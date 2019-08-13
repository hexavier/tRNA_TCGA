library(ggplot2)
library(ggpubr)

# Load nTE data
load(file="/home/xhernandez/Documents/tAI_genomic/AAcTE_CUgenomic_sqrt.rda")

# Rename samples based on clustering
colnames(nTE) = gsub("KICH","kidney",colnames(nTE))
colnames(nTE) = gsub("KIRP","kidney",colnames(nTE))
colnames(nTE) = gsub("KIRC","kidney",colnames(nTE))
colnames(nTE) = gsub("LUAD","lung",colnames(nTE))
colnames(nTE) = gsub("LUSC","lung",colnames(nTE))
colnames(nTE) = gsub("COAD","colorectal",colnames(nTE))
colnames(nTE) = gsub("READ","colorectal",colnames(nTE))
colnames(nTE) = gsub("CHOL","liver",colnames(nTE))
colnames(nTE) = gsub("LIHC","liver",colnames(nTE))
colnames(nTE) = gsub("UCEC","uterus",colnames(nTE))
colnames(nTE) = gsub("CESC","uterus",colnames(nTE))

rownames(nTE) = gsub("KICH","kidney",rownames(nTE))
rownames(nTE) = gsub("KIRP","kidney",rownames(nTE))
rownames(nTE) = gsub("KIRC","kidney",rownames(nTE))
rownames(nTE) = gsub("LUAD","lung",rownames(nTE))
rownames(nTE) = gsub("LUSC","lung",rownames(nTE))
rownames(nTE) = gsub("COAD","colorectal",rownames(nTE))
rownames(nTE) = gsub("READ","colorectal",rownames(nTE))
rownames(nTE) = gsub("CHOL","liver",rownames(nTE))
rownames(nTE) = gsub("LIHC","liver",rownames(nTE))
rownames(nTE) = gsub("UCEC","uterus",rownames(nTE))
rownames(nTE) = gsub("CESC","uterus",rownames(nTE))

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
pdf("plots/AAcTE_bycod_codons_unpaired.pdf",height=5,width = 5)
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
