library(gplots)
library(ggplot2)

#%% Load data
cancer_types = c(BRCA="BRCA",PRAD="PRAD",kidney="KICH;KIRP;KIRC",lung="LUAD;LUSC",HNSC="HNSC",uterus="UCEC;CESC",
                 liver="LIHC;CHOL",THCA="THCA",colorectal="COAD;READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG")
cancer_types = c(BRCA="BRCA",PRAD="PRAD",KICH="KICH",KIRC="KIRC",KIRP="KIRP",LUAD="LUAD",LUSC="LUSC",HNSC="HNSC",UCEC="UCEC",CESC="CESC",
                 LIHC="LIHC", CHOL="CHOL",THCA="THCA",COAD="COAD",READ="READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG")
anticodon = read.csv("results/AAcTE_CUgenomic_sqrt.csv",row.names = 1)
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)
rownames(anticodon) = paste0(codons[rownames(anticodon),"AA"],rownames(anticodon))

# Initialize structures
diffexp_mean = matrix(nrow = nrow(anticodon),ncol = length(cancer_types))
rownames(diffexp_mean) = rownames(anticodon); colnames(diffexp_mean) = names(cancer_types)
is_sig = matrix(nrow = nrow(anticodon),ncol = length(cancer_types))
rownames(is_sig) = rownames(anticodon); colnames(is_sig) = names(cancer_types)

# Scatter plot deltaPSI vs SD
for (type in names(cancer_types)){
  idxtissue = (sapply(colnames(anticodon), function(x) strsplit(x,"\\.")[[1]][1]) %in% strsplit(cancer_types[type],";")[[1]])
  # Upload gene expression and methylation data
  tempdata = anticodon[,idxtissue]
  colnames(tempdata) = sapply(colnames(tempdata), function(x) paste(strsplit(x,"\\.")[[1]][2:10], collapse = "-"))
  
  ## Analyze methylation status of genes
  expr_samples = substr(colnames(tempdata),1,15)
  is_normal = (regexpr("TCGA-[A-Z,0-9]{2}-[A-Z,0-9]{4}-11",expr_samples))==1
  is_tum = (regexpr("TCGA-[A-Z,0-9]{2}-[A-Z,0-9]{4}-01",expr_samples))==1
  
  # Add data
  diffexp_mean[,type] = log2(as.numeric(rowMeans(tempdata[,is_tum],na.rm = T)) / 
                                 as.numeric(rowMeans(tempdata[,is_normal],na.rm = T)))
  # Significance
  pvals = sapply(rownames(tempdata), function(x) wilcox.test(as.numeric(tempdata[x,is_tum]),
                                        as.numeric(tempdata[x,is_normal]))$p.value)

  is_sig[,type] = pvals
}

# Multiple testing correction
is_sig[,] <- p.adjust(is_sig,method="fdr")

# Keep only significant exons
diffexp_mean[(is_sig>=0.05)|(is.na(is_sig))] = NA

# Heatmap
breaks = c(min(diffexp_mean,na.rm=T),seq(-0.5,0.5,length.out=254),max(diffexp_mean,na.rm=T))
diffexp_mean = diffexp_mean[order((rowSums(diffexp_mean<0,na.rm=T)-rowSums(diffexp_mean>0,na.rm=T)),decreasing = T),order(colSums(is.na(diffexp_mean)))]
heatmap.2(diffexp_mean,
          Rowv=F, Colv=F, col = bluered(255), symm=T, margins=c(9,15), trace="none",
          ylab="delta PSI", xlab="Cancer Types", na.color="grey",breaks = breaks)

# Heatmap based on healthy pca weights
AAcTE_pca = read.csv("results/AAcTE_pca.txt", row.names = 1, sep = "\t")
breaks = c(min(diffexp_mean,na.rm=T),seq(-0.5,0.5,length.out=254),max(diffexp_mean,na.rm=T))
diffexp_mean=diffexp_mean[order(AAcTE_pca[rownames(diffexp_mean),1],decreasing = T),order(colSums(is.na(diffexp_mean)))]
heatmap.2(diffexp_mean,
          Rowv=F, Colv=F, col = bluered(255), symm=T, margins=c(9,15), trace="none",
          ylab="delta PSI", xlab="Cancer Types", na.color="grey",breaks = breaks)

# Two sided bar plot to show # cancer types in each direction
cancer_counts = data.frame(row.names=1:(nrow(diffexp_mean)*2))
cancer_counts$num = c(rowSums(diffexp_mean>0, na.rm = T),- rowSums(diffexp_mean<0, na.rm = T))
cancer_counts$group = c(rep("UP", nrow(diffexp_mean)), rep("DOWN", nrow(diffexp_mean)))
cancer_counts$order = rep(1:(nrow(diffexp_mean)),2)

ggplot(cancer_counts, aes(x = order, y = num, fill = group)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("blue","red")) + 
  theme_minimal()