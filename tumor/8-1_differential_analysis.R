library(gplots)

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

# Build summary table statistics
summary_table = data.frame(row.names = 1:(length(cancer_types)*nrow(anticodon)))
catypes = c()
isoacceptor = c()
meansT = c()
meansH = c()
pvals = c()
for (g in rownames(anticodon)){
  catypes = append(catypes,names(cancer_types))
  isoacceptor = append(isoacceptor, rep(g,length(cancer_types)))
  meansT = append(meansT,sapply(names(cancer_types),
                                    function(x,dataset) mean(dataset[(dataset$catype==x)&(dataset$state=="CA"),g],na.rm=T),dataset))
  meansH = append(meansH,sapply(names(cancer_types),
                                    function(x,dataset) mean(dataset[(dataset$catype==x)&(dataset$state=="HE"),g],na.rm=T),dataset))
  pvals = append(pvals,sapply(names(cancer_types),
                              function(x,dataset) wilcox.test(dataset[(dataset$catype==x)&(dataset$state=="HE"),g],
                                                              dataset[(dataset$catype==x)&(dataset$state=="CA"),g],
                                                              alternative = "two.sided")$p.value,dataset))
}
summary_table$catype = catypes
summary_table$isoacceptor = isoacceptor
summary_table$mean_exp_T = meansT
summary_table$mean_exp_H = meansH
summary_table$FC = log2(meansT/meansH)
summary_table$pval = pvals
summary_table$pval_adj = p.adjust(pvals,method="fdr")
write.csv(summary_table, sprintf("results/differential_AAcTE_TvsH_twosided.csv",type))