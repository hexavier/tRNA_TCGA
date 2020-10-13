library(gplots)

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
  ## Analyze methylation status of genes
  mapped_normal = (regexpr(sprintf("TCGA_%s_N{1}[A-Z,0-9]{4}",type),colnames(rawvalues)))==1
  
  dataset_temp = data.frame(row.names = colnames(rawvalues))
  dataset_temp$sample = colnames(rawvalues)
  dataset_temp[,rownames(rawvalues)] = t(rawvalues)
  dataset_temp$catype = type
  dataset_temp$state = "CA"
  dataset_temp[mapped_normal,"state"] = "HE"
  dataset = rbind(dataset, dataset_temp)
}

# Build summary table statistics
cancer_types = unique(types)

summary_table = data.frame(row.names = 1:(length(cancer_types)*nrow(rawvalues)))
catypes = c()
gene = c()
isoacceptor = c()
meansT = c()
meansH = c()
pvals = c()
for (g in rownames(rawvalues)){
  catypes = append(catypes,cancer_types)
  gene = append(gene, rep(g,length(cancer_types)))
  isoacceptor = append(isoacceptor, rep(strsplit(g,"-")[[1]][2],length(cancer_types)))
  meansT = append(meansT,sapply(cancer_types,
                                function(x,dataset) mean(dataset[(dataset$catype==x)&(dataset$state=="CA"),g],na.rm=T),dataset))
  meansH = append(meansH,sapply(cancer_types,
                                function(x,dataset) mean(dataset[(dataset$catype==x)&(dataset$state=="HE"),g],na.rm=T),dataset))
}
summary_table$catype = catypes
summary_table$gene = gene
summary_table$isoacceptor = isoacceptor
summary_table$mean_exp_T = meansT
summary_table$mean_exp_H = meansH
summary_table$delta = (meansT-meansH)
write.csv(summary_table, sprintf("results/differentialBisulfite_TvsH_twosided.csv",type))

#### CONSISTENCY ####
values=summary_table
#Initialize structure
anticodons = unique(values$isoacceptor)
consistency = matrix(nrow = length(anticodons), ncol=length(cancer_types))
rownames(consistency) = anticodons; colnames(consistency) = cancer_types
for (c in anticodons){
  aadata = values[values$isoacceptor==c,]
  up = sapply(cancer_types,function(x) sum(aadata[aadata$catype==x,"delta"]>0,na.rm = T)/sum(!is.na(aadata[aadata$catype==x,"delta"])))
  down = sapply(cancer_types,function(x) sum(aadata[aadata$catype==x,"delta"]<0,na.rm = T)/sum(!is.na(aadata[aadata$catype==x,"delta"])))
  consistency[c,] = (up-down)
}

# Heatmap
breaks = c(min(consistency,na.rm=T),seq(-0.5,0.5,length.out=254),max(consistency,na.rm=T))
diffexp_mean = consistency[order((rowSums(consistency<0,na.rm=T)-rowSums(consistency>0,na.rm=T)),decreasing = T),order(colSums(is.na(consistency)))]
heatmap.2(diffexp_mean,
          Rowv=F, Colv=F, col = bluered(255), symm=T, margins=c(12,15), trace="none",
          na.color="grey",breaks = breaks)

#### HEATMAP INTERESTING GENES ####
aadata = values[substr(values$isoacceptor,1,3)=="Arg",]
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
breaks = c(min(diffexp_mean,na.rm=T),seq(-30,30,length.out=254),max(diffexp_mean,na.rm=T))
diffexp_mean = diffexp_mean[order((rowSums(diffexp_mean<0,na.rm=T)-rowSums(diffexp_mean>0,na.rm=T)),decreasing = T),order(colSums(is.na(diffexp_mean)))]
heatmap.2(diffexp_mean,
          Rowv=F, Colv=F, col = bluered(255), symm=T, margins=c(12,15), trace="none",
          na.color="grey",breaks = breaks)