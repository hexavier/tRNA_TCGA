library(ggplot2)
library(ggpubr)

extract_cod <- function (trnas, anticod){
  output = data.frame(row.names = anticod)
  trnas_acod = sapply(rownames(trnas), function(x) substr(x,nchar(x)-2,nchar(x)))
  for (s in colnames(trnas)){
    output[,s] = sapply(anticod, function(x) sum(trnas[trnas_acod %in% x,s]))
  }
  return(output)
}

transformdata <- function(data,transf){
  aa_idx = regexpr("i?[A-Z][a-z]{2}[A-Z]{3}",rownames(data))==1
  data = data[aa_idx,]
  if (transf=="log"){
    outdata = sapply(data,log)
    # Remove inf values
    outdata[outdata==-Inf] = NaN
    rownames(outdata)=rownames(data)
  }else if (transf=="arcsinh"){
    outdata = sapply(data,asinh)
    rownames(outdata)=rownames(data)
  }else if (transf=="sqrt"){
    outdata = sapply(data,sqrt)
    rownames(outdata)=rownames(data)
  }else if (transf=="rel"){
    # Compute relative data
    outdata = data.frame(matrix(ncol = ncol(data), nrow = nrow(data)),row.names = rownames(data))
    colnames(outdata)= colnames(data)
    aa = sapply(rownames(outdata),function(x) substr(x,1,nchar(x)-3))
    uniqueaa = unique(aa)
    for (n in uniqueaa){
      idx = (aa %in% n)
      idx_data = matrix(as.matrix(data[idx,]), ncol = ncol(data), nrow = sum(idx))
      total = colSums(idx_data)
      outdata[idx,] = t(apply(idx_data,1,function(x) x/total))
      iszero = (total %in% 0)
      if (any(iszero)){
        outdata[idx,iszero] = 1.0/sum(idx)
      }
    }
  }else{
    outdata=data
  }
  return(outdata)
}

# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# tRNAs
trna = read.csv("data/HannahSeq_nomod.csv",row.names = 1)
anticodon = extract_cod(transformdata(trna,"sqrt"),codons$ANTICODON)

# See cluster of main differences
clusters = hclust(dist(t(anticodon))) # samples
plot(clusters)

# To see correlations
cortable = data.frame(cor(anticodon, method="spearman")[1:15,16:ncol(anticodon)])
labels = as.character(sapply(rownames(cortable),function(x) strsplit(x,"\\.")[[1]][1]))

# Reestructure data in rows
dataset = c()

for (sample in colnames(cortable)){
  dataset_temp = data.frame(row.names=seq(nrow(cortable)))
  dataset_temp$hydro = strsplit(sample,"\\.")[[1]][1]
  dataset_temp$correlation = as.numeric(cortable[,sample])
  dataset_temp$match = (labels==(strsplit(sample,"\\.")[[1]][1]))
  dataset = rbind(dataset, dataset_temp)
}


# Differential expression
diffexp = compare_means(correlation~match, data=dataset, group.by = "hydro", method="wilcox.test",alternative = "greater")
diffexp$log2FC = apply(diffexp,1, function(x) log2(mean(dataset[(dataset$match==x["group1"])&(dataset$hydro==x["hydro"]),"correlation"])/
                                                 mean(dataset[(dataset$match==x["group2"])&(dataset$hydro==x["hydro"]),"correlation"])))

# Plot
ggplot(dataset, aes(x=hydro, y=correlation, fill=match)) +  
  geom_violin(position=position_dodge(1)) +
  labs(title="Hydro-seq VS small RNA-seq",x="Hydro-seq sample", y = "Spearman correlation") + 
  stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
  stat_compare_means(aes(label = ..p.format..), method="wilcox.test",method.args = list(alternative = "greater")) +
  theme_classic()
