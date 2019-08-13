library(gplots)

# Load nTE data
load(file="/home/xhernandez/Documents/tAI_genomic/AAcTE_CUgenomic_sqrt.rda")
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# Extract matched samples
matched = apply(nTE,3,diag)
colnames(matched) = paste0(codons[colnames(matched),"AA"],colnames(matched))
write.csv(matched,"results/matched_AAcTE.csv")

# Restructure data
dataset = data.frame()
tissues = sapply(rownames(matched), function(x) strsplit(x,"\\.")[[1]][1])
for (cod in colnames(matched)){
  data_temp = data.frame(row.names = rownames(matched))
  data_temp$nTE = matched[,cod]
  data_temp$codon = cod
  data_temp$tissue = tissues
  dataset = rbind(dataset,data_temp)
}
colnames(dataset) = c("nTE","codon","tissue")

# ANOVA
# nTE_aov = aov(nTE ~ codon*tissue, data=dataset)
# summary(nTE_aov)
# plot(nTE_aov)

## Plot
medians = matrix(nrow=ncol(matched),ncol=length(unique(tissues)))
rownames(medians)=colnames(matched); colnames(medians) = unique(tissues)
for (t in colnames(medians)){
  medians[,t]=apply(matched[tissues==t,],2,median,na.rm=T)
}
medians = medians[order(rowMeans(medians)),]

heatmap.2(medians,symm=F)