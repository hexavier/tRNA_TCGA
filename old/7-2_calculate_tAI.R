library(tAI)

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

## Load trna and weighted CU
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# tRNAs
trna = read.csv("data/TCGAall_nomod.csv",row.names = 1)
anticodon = extract_cod(transformdata(trna,"arcsinh"),codons$ANTICODON)

# Codon usage
weighted_CU = read.csv("results/healthy_weightedCU.csv",row.names = 1)
short_samples = sapply(colnames(trna), function(x) paste(strsplit(x,"\\.")[[1]][2:5],collapse="."))
matched_CU = data.frame(sapply(short_samples, function(x) if(any((substr(colnames(weighted_CU),1,16) %in% x))){
  weighted_CU[,(substr(colnames(weighted_CU),1,16) %in% x)]
  }else{matrix(data=NA,nrow=1,ncol=64)}))
rownames(matched_CU) = sapply(rownames(weighted_CU),function(x) paste(codons[x,"AA"],x,sep=""))
codon = extract_cod(transformdata(matched_CU,""),rownames(codons)[!(codons$AA %in% c("Stop","Met"))])

TAIs = data.frame(matrix(ncol = ncol(anticodon), nrow = ncol(anticodon)),row.names = colnames(anticodon)); colnames(TAIs) = colnames(anticodon)
initial_s = c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)
## Calculate tAI
for (sample in colnames(trna)){
  # Calculate relative adaptiveness values (ws)
  sample.ws = get.ws(tRNA=anticodon[,sample], s=initial_s, sking=0)
  # Calculate tAI for all CUs
  sample.tai <- get.tai(t(codon), sample.ws)
  TAIs[sample,] = sample.tai
}

write.csv(TAIs,"results/tAI_weightedcodons.csv")