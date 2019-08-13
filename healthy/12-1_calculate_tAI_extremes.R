library(tAI)
library(ggplot2)

extract_cod <- function (trnas, anticod){
  output = data.frame(row.names = anticod)
  trnas_acod = sapply(rownames(trnas), function(x) substr(x,nchar(x)-2,nchar(x)))
  for (s in colnames(trnas)){
    output[,s] = sapply(anticod, function(x) if(any(trnas_acod==x)){mean(trnas[trnas_acod==x,s])}else{0})
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

AAnormalize <- function(data,codons){
  # Compute relative data
  aa = sapply(codons,function(x) substr(x,1,nchar(x)-3))
  uniqueaa = unique(aa)
  outdata = numeric(length=length(data))
  for (n in uniqueaa){
    idx = (aa %in% n)
    total = max(data[idx],na.rm=T)
    outdata[idx] = data[idx]/total
    if (total %in% 0){
      outdata[idx] = 1.0/sum(idx)
    }
  }
  return(outdata)
}

## Load trna and weighted CU
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# tRNAs
trna = read.csv("data/TCGAall_nomod.csv",row.names = 1)
colnames(trna) = gsub("COAD","colorectal",colnames(trna))
colnames(trna) = gsub("READ","colorectal",colnames(trna))
# Keep only extremes
tissues = sapply(colnames(trna), function(x) strsplit(x, "\\.")[[1]][1])
trna = trna[,(tissues %in% c("GBM","colorectal"))]
# Compute mean of tissue
tissues = sapply(colnames(trna), function(x) strsplit(x, "\\.")[[1]][1])
trna_mean = data.frame(row.names = rownames(trna))
trna_mean$GBM = rowMeans(trna[,tissues=="GBM"])
trna_mean$colorectal = rowMeans(trna[,tissues=="colorectal"])
anticodon = extract_cod(transformdata(trna_mean,"sqrt"),codons$ANTICODON)

# Genomic codon usage
codus = read.csv("data/human_CU_refseq_CoCoPUT.tsv",sep="\t")
codus_idx = as.character(codus$Protein.ID)
# Keep only columns with codon info
codus_clean = data.frame(sapply(unique(codus_idx),function(x) colMeans(codus[(codus_idx==x),13:ncol(codus)],na.rm = T)))
rownames(codus_clean) = sapply(rownames(codus_clean),function(x) paste(codons[x,"AA"],x,sep=""))

## Calculate tAI for genomic CU
codon = extract_cod(transformdata(codus_clean,""),rownames(codons)[!(codons$AA %in% c("Stop","Met"))])

TAIs = data.frame(matrix(ncol = ncol(anticodon), nrow = ncol(codon)),row.names = colnames(codon)); colnames(TAIs) = colnames(anticodon)
initial_s = c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)
# Calculate tAI
for (sample in colnames(anticodon)){
  # Calculate relative adaptiveness values (ws)
  sample.ws = get.ws(tRNA=anticodon[,sample], s=initial_s, sking=0)
  #sample.ws = AAnormalize(sample.ws,paste0(codons[rownames(codon),"AA"],rownames(codon))) # normalize by AA
  # Calculate tAI for all CUs
  sample.tai <- get.tai(t(codon), sample.ws)
  TAIs[,sample] = sample.tai
}
write.csv(TAIs,"results/tAI_2extremes.csv")