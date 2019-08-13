library(gplots)

#%% Load data
anticodon = read.csv("results/AAcTE_CUproteomics_sqrt.csv",row.names = 1)
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)
rownames(anticodon) = paste0(codons[rownames(anticodon),"AA"],rownames(anticodon))

## Plot
tissues = sapply(colnames(anticodon), function(x) strsplit(x,"\\.")[[1]][1])
medians = matrix(nrow=nrow(anticodon),ncol=length(unique(tissues)))
rownames(medians)=rownames(anticodon); colnames(medians) = unique(tissues)
for (t in colnames(medians)){
  medians[,t]=apply(anticodon[,tissues==t],1,median,na.rm=T)
}
medians = medians[order(rowMeans(medians)),]

heatmap.2(medians,symm=F,cexCol = 1)
