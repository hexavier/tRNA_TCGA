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
    #total = sum(data[idx],na.rm=T)
    total = max(data[idx],na.rm=T)
    outdata[idx] = data[idx]/total
    if (total %in% 0){
      outdata[idx] = 1.0/sum(idx)
    }
  }
  return(outdata)
}

scaling01 <- function(df){
  df_out = apply(df,2,function(x) (x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T)))
  return(df_out)
}

# Define cancer types
cancer_types = c(BRCA="BRCA",colorectal="COAD;READ")

## Load trna and weighted CU
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# tRNAs
trnaH = read.csv("data/TCGAhealthy_nomod.csv",row.names = 1)
trnaT = read.csv("data/TCGAtumor_nomod.csv",row.names = 1)
trna = cbind(trnaH,trnaT)
anticodon = extract_cod(transformdata(trna,"sqrt"),codons$ANTICODON)

# Genomic codon usage
codus = read.csv("data/human_CU_refseq.tsv",sep="\t")

# Mapping
np2nm = read.csv("data/NP2NM.txt",sep="\t", header=T)
#http://genome-euro.ucsc.edu/cgi-bin/hgTables ncbiRefLink
np2nm = np2nm[np2nm[,6]!="",]
rownames(np2nm) = sapply(np2nm[,6], function(x) substr(as.character(x),1,nchar(as.character(x))-2))

initial_s = c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)
for (type in names(cancer_types)){
  # Protein gene expression
  protseq = read.csv(sprintf("/home/xhernandez/Downloads/TCGA-proteomics/TCGA_%s.csv",type),row.names = 1)
  idxtissue = (sapply(colnames(anticodon), function(x) strsplit(x,"\\.")[[1]][1]) %in% strsplit(cancer_types[type],";")[[1]])
  # Remove other tissues
  temptrna = anticodon[,idxtissue]
  # Map gene expression samples and remove others
  short_samples = sapply(colnames(temptrna), function(x) paste(strsplit(x,"\\.")[[1]][2:5],collapse="."))
  short_genexp_samples = substr(colnames(protseq),1,16)
  merged_samples = short_samples[short_samples %in% short_genexp_samples]
  temptrna = temptrna[,short_samples %in% short_genexp_samples]
  
  # Calculate normalized TE based on demand and supply of codons
  nTE = matrix(NaN,ncol = ncol(temptrna), nrow = length(rownames(codons)[!(codons$AA %in% c("Stop","Met"))]));
  colnames(nTE) = colnames(temptrna); rownames(nTE) = rownames(codons)[!(codons$AA %in% c("Stop","Met"))]
  
  ## Prepare CU data
  # Map ids
  gene_ids = rownames(protseq)[(rownames(protseq) %in% np2nm$name)]
  refseq_ids = sapply(gene_ids, function(x) sapply(np2nm[(np2nm$name %in% x),"protAcc"],function(s)substr(s,1,nchar(as.character(s))-2)))
  # Keep only protAcc that qre present in codus
  codus_ids = sapply(codus$Protein.ID,function(x) substr(x,1,nchar(as.character(x))-2))
  refseq_ids = lapply(refseq_ids, function(x) x[x %in% codus_ids])
  # Convert booleans in codus_ids to string index. If more than 1, take mean. Remove NA values(not present in codus in last step)
  codus_clean = data.frame(lapply(refseq_ids,function(x) colMeans(codus[(codus_ids %in% x),13:ncol(codus)],na.rm = T)))
  codus_clean = codus_clean[,!is.na(codus_clean[1,])]
  
  # Prepare codon data
  rownames(codus_clean) = sapply(rownames(codus_clean),function(x) paste(codons[x,"AA"],x,sep=""))
  codon = extract_cod(transformdata(codus_clean,""),rownames(codons)[!(codons$AA %in% c("Stop","Met"))])
  
  # Clean isofseq dataset
  protseq = protseq[colnames(codon),]
  # Transform values to 0-1
  protseq = scaling01(protseq)
  ################################################ RANDOMIZE ######################################################
  protseq = apply(protseq,2,sample)
  #################################################################################################################
  
  # Calculate expression of matched samples
  matched_genexp = data.frame(sapply(merged_samples, function(x) if(sum(short_genexp_samples==x)==1){protseq[,(short_genexp_samples==x)]
  }else{rowMeans(protseq[,(short_genexp_samples==x)],na.rm = T)}))
  # Do matrix multiplication, substituting NA to 0 in order to skip non detected proteins
  matched_genexp[is.na(matched_genexp)] <- 0
  CU = apply(matched_genexp,2,function(x,codon) as.matrix(codon)%*%as.matrix(x),codon)
  CU = t(apply(CU,2,AAnormalize,paste0(codons[rownames(codon),"AA"],rownames(codon)))) # normalize by AA
  
  # Calculate tAI
  sample.ws = apply(temptrna,2, get.ws, s=initial_s, sking=0)
  sample.ws = t(apply(sample.ws,2,AAnormalize,paste0(codons[rownames(codon),"AA"],rownames(codon)))) # normalize by AA
  
  # Compute ratio between supply and demand
  cTE = sample.ws/CU
  nTE[,] = t(cTE) # keep non-normalized results
  write.csv(nTE,sprintf("results/AAcTE_RANDOM_CUproteomics%s_sqrt.csv",type))
}