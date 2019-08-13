library(ggplot2)
library(ggpubr)

#%% Load data
cancer_types = c(BRCA="BRCA",PRAD="PRAD",KICH="KICH",KIRC="KIRC",KIRP="KIRP",LUAD="LUAD",LUSC="LUSC",HNSC="HNSC",UCEC="UCEC",CESC="CESC",
                 LIHC="LIHC", CHOL="CHOL",THCA="THCA",COAD="COAD",READ="READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG")

anticodon = read.csv("results/AAcTE_CUgenomic_sqrt.csv",row.names = 1)

# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)
rownames(anticodon) = paste0(codons[rownames(anticodon),"AA"],rownames(anticodon))

## Differential analysis
dataset = c()

for (type in names(cancer_types)){
  idxtissue = (sapply(colnames(anticodon), function(x) strsplit(x,"\\.")[[1]][1]) %in% strsplit(cancer_types[type],";")[[1]])
  # Upload gene expression and methylation data
  tempdata = anticodon[,idxtissue]
  colnames(tempdata) = sapply(colnames(tempdata), function(x) paste(strsplit(x,"\\.")[[1]][2:10], collapse = "-"))
  
  ## Analyze methylation status of genes
  expr_samples = substr(colnames(tempdata),1,15)
  mapped_normal = (regexpr("TCGA-[A-Z,0-9]{2}-[A-Z,0-9]{4}-11",expr_samples))==1
  
  dataset_temp = data.frame(row.names = colnames(tempdata))
  dataset_temp$sample = expr_samples
  dataset_temp[,rownames(tempdata)] = t(tempdata)
  dataset_temp$catype = type
  dataset_temp$state = "CA"
  dataset_temp[mapped_normal,"state"] = "HE"
  dataset = rbind(dataset, dataset_temp)
}

pdf("plots/AAcTE_diffexp_boxplots_rel.pdf",height=5, width=12)
for (gene in rownames(anticodon)){
  print(ggplot(dataset, aes(x=catype, y=get(gene), fill=state)) +  
          geom_boxplot(position=position_dodge(1)) +
          scale_fill_manual(values=c("red", "blue")) +
          labs(title=sprintf("%s expression across TCGA",gene),x="Cancer Type", y = "Gene Expression") + 
          stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
          stat_compare_means(aes(label = ..p.signif..)) +
          theme_classic())
}
dev.off()