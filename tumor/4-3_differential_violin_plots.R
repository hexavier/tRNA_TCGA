library(ggplot2)
library(ggpubr)

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

#%% Load data
cancer_types = c(BRCA="BRCA",PRAD="PRAD",KICH="KICH",KIRC="KIRC",KIRP="KIRP",LUAD="LUAD",LUSC="LUSC",HNSC="HNSC",UCEC="UCEC",CESC="CESC",
                 LIHC="LIHC", CHOL="CHOL",THCA="THCA",COAD="COAD",READ="READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG")

trnaH = read.csv("data/TCGAhealthy_nomod.csv",row.names = 1)
trnaT = read.csv("data/TCGAtumor_nomod.csv",row.names = 1)
trna = cbind(trnaH,trnaT)
anticodon = transformdata(trna,"")

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

pdf("plots/diffexp_boxplots.pdf",height=8, width=12)
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