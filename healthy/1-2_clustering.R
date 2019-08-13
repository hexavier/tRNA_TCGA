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
anticodon = transformdata(trna,"sqrt")

# Calculate median per tissue
tissues = sapply(colnames(anticodon), function(x) strsplit(x,"\\.")[[1]][1])
medians = matrix(nrow=nrow(anticodon),ncol=length(unique(tissues)))
rownames(medians)=rownames(anticodon); colnames(medians) = unique(tissues)
for (t in colnames(medians)){
  medians[,t]=apply(anticodon[,tissues==t],1,median,na.rm=T)
}
medians = medians[order(rowMeans(medians)),]

# Heatmap
## Plot
heatmap.2(medians,symm=F)

