library(gplots)
library(ggplot2)

mean_anticodon <- function(df){
  anticodons = sapply(rownames(df),function(x) strsplit(x,"-")[[1]][2])
  df_out = t(sapply(unique(anticodons), function(x) colMeans(df[anticodons==x,], na.rm=T)))
  return(df_out)
}

# Load data
path="/users/lserrano/xhernandez/tRNA_methylation/"

## Differential analysis Methylation
allvalues = read.csv(sprintf("%sResults/trnas_bisulfite.csv",path), row.names = 1)
# Calculate % coverage
percent = apply(allvalues,2,function(x) sum(!is.na(x))/length(x))

dataset = c()

types = sapply(colnames(allvalues), function(x) strsplit(x,"_")[[1]][2])
for (type in unique(types)){
  # Get type samples
  rawvalues = allvalues[,(types %in% type)]
  meanvalues = mean_anticodon(rawvalues)
  ## Analyze methylation status of genes
  mapped_normal = (regexpr(sprintf("TCGA_%s_N{1}[A-Z,0-9]{4}",type),colnames(meanvalues)))==1
  
  dataset_temp = data.frame(row.names = colnames(meanvalues))
  dataset_temp$sample = colnames(meanvalues)
  dataset_temp[,rownames(meanvalues)] = t(meanvalues)
  dataset_temp$catype = type
  dataset_temp$state = "CA"
  dataset_temp[mapped_normal,"state"] = "HE"
  dataset = rbind(dataset, dataset_temp)
}

# Build summary table statistics
cancer_types = unique(types)

summary_table = data.frame(row.names = 1:(length(cancer_types)*nrow(meanvalues)))
catypes = c()
gene = c()
meansT = c()
meansH = c()
pvals = c()
for (g in rownames(meanvalues)){
  catypes = append(catypes,cancer_types)
  gene = append(gene, rep(g,length(cancer_types)))
  meansT = append(meansT,sapply(cancer_types,
                                function(x,dataset) median(dataset[(dataset$catype==x)&(dataset$state=="CA"),g],na.rm=T),dataset))
  meansH = append(meansH,sapply(cancer_types,
                                function(x,dataset) median(dataset[(dataset$catype==x)&(dataset$state=="HE"),g],na.rm=T),dataset))
}
summary_table$catype = catypes
summary_table$gene = gene
summary_table$mean_exp_T = meansT
summary_table$mean_exp_H = meansH
summary_table$delta = (meansT-meansH)
write.csv(summary_table, sprintf("results/differentialBisulfite_MEDIANS_TvsH_twosided.csv",type))


#### HEATMAP INTERESTING GENES ####
values=summary_table
#aadata = values[substr(values$gene,1,3)=="Pro",]
aadata=values
genes = unique(aadata$gene)

# Initialize structures
diffexp_mean = matrix(nrow = length(genes), ncol = length(cancer_types))
rownames(diffexp_mean) = genes; colnames(diffexp_mean) = cancer_types

for (type in cancer_types){
  # Add data
  diffexp_mean[,type] = sapply(genes, function(x) aadata[(aadata$catype==type)&(aadata$gene==x),"delta"])
}

# Keep only significant exons
diffexp_mean = diffexp_mean[rowSums(!is.na(diffexp_mean))>0,]

# Heatmap
breaks = c(min(diffexp_mean,na.rm=T),seq(-20,20,length.out=254),max(diffexp_mean,na.rm=T))
diffexp_mean = diffexp_mean[order((rowSums(diffexp_mean<0,na.rm=T)-rowSums(diffexp_mean>0,na.rm=T)),decreasing = T),order(colSums(is.na(diffexp_mean)))]
heatmap.2(diffexp_mean,
          Rowv=F, Colv=F, col = bluered(255), symm=T, margins=c(8,8), trace="none",
          na.color="grey",breaks = breaks)

# Two sided bar plot to show # cancer types in each direction
cancer_counts = data.frame(row.names=1:(nrow(diffexp_mean)*2))
cancer_counts$num = c(rowSums(diffexp_mean>0, na.rm = T),- rowSums(diffexp_mean<0, na.rm = T))
cancer_counts$group = c(rep("UP", nrow(diffexp_mean)), rep("DOWN", nrow(diffexp_mean)))
cancer_counts$order = rep(1:(nrow(diffexp_mean)),2)

ggplot(cancer_counts, aes(x = order, y = num, fill = group)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("blue","red")) + 
  theme_minimal()