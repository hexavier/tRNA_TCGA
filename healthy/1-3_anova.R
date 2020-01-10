library(gplots)

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

# tRNAs
trna = read.csv("data/TCGAall_nomod.csv",row.names = 1)
anticodon = transformdata(trna,"")

# Re-structure dataset
tissues = sapply(colnames(anticodon), function(x) strsplit(x,"\\.")[[1]][1])
dataset = c()
for (s in colnames(anticodon)){
  dataset_temp = data.frame(row.names = 1:nrow(anticodon))
  dataset_temp$expr = as.numeric(anticodon[,s])
  dataset_temp$codon = rownames(anticodon)
  dataset_temp$tissue = strsplit(s,"\\.")[[1]][1]
  dataset = rbind(dataset,dataset_temp)
}

# Data does not match assumptions
res.kruskal = list()
for (c in unique(dataset$codon)){
  res.kruskal[[c]] = kruskal.test(expr ~ tissue, data=dataset[dataset$codon==c,])
}
pvals = as.numeric(lapply(res.kruskal,function(x) unlist(x[[3]])))
corr_pvals = p.adjust(pvals,method = "fdr")
