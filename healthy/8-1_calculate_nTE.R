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

## Load trna and weighted CU
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# tRNAs
trna = read.csv("data/TCGAall_nomod.csv",row.names = 1)
anticodon = extract_cod(transformdata(trna,"sqrt"),codons$ANTICODON)

# Genomic codon usage
codus = read.csv("data/human_CU_refseq_CoCoPUT.tsv",sep="\t")

# Isoform gene expression
isofseq = read.csv("data/healthy_isoform_seq.csv",row.names = 1)
rownames(isofseq) = as.character(sapply(rownames(isofseq), function(x) substr(x,1,nchar(x)-2)))

# Mapping
tcga2refseq = read.csv("data/TCGA2RefSeq.txt",sep="\t", header=F,row.names = 1)
# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownIsoforms.txt.gz
rownames(tcga2refseq) = sapply(rownames(tcga2refseq), function(x) substr(x,1,nchar(x)-2))

np2nm = read.csv("data/NP2NM.txt",sep="\t", header=T)
#http://genome-euro.ucsc.edu/cgi-bin/hgTables ncbiRefLink
np2nm = np2nm[np2nm[,6]!="",]
rownames(np2nm) = sapply(np2nm[,6], function(x) substr(as.character(x),1,nchar(as.character(x))-2))

## Prepare CU data
# Map ids
tcga_ids = rownames(isofseq)[(rownames(isofseq) %in% rownames(tcga2refseq))]
refseq_ids = sapply(rownames(isofseq), function(x) if (x %in% rownames(tcga2refseq)){as.character(tcga2refseq[x,1])}else{NA})
refseq_ids = as.character(refseq_ids[!is.na(refseq_ids)])
# Detect whether refseq ids are in codon usage data, and remove if not
codusindex_main = sapply(codus$Protein.ID, function(x) substr(as.character(np2nm[substr(x,1,nchar(as.character(x))-2),"X.id"]),1,nchar(as.character(x))-2))
# Keep only events in common between tcga and codon usage
keep = (refseq_ids %in% codusindex_main)
tcga_ids = tcga_ids[keep]
refseq_ids = refseq_ids[keep]
# Convert booleans in codus_ids to string index. If more than 1, take mean
codus_clean = data.frame(sapply(refseq_ids,function(x) colMeans(codus[(codusindex_main==x),13:ncol(codus)],na.rm = T)))

# Clean isofseq dataset
isofseq = isofseq[tcga_ids,]

# Prepare codon data
rownames(codus_clean) = sapply(rownames(codus_clean),function(x) paste(codons[x,"AA"],x,sep=""))
codon = extract_cod(transformdata(codus_clean,""),rownames(codons)[!(codons$AA %in% c("Stop","Met"))])

# Map gene expression samples
short_samples = sapply(colnames(anticodon), function(x) paste(strsplit(x,"\\.")[[1]][2:5],collapse="."))
short_genexp_samples = substr(colnames(isofseq),1,16)
matched_genexp = data.frame(sapply(short_samples, function(x) if(any((short_genexp_samples==x))){
  isofseq[,(short_genexp_samples==x)]}else{matrix(data=NA,nrow=1,ncol=ncol(codon))}))
CU = apply(matched_genexp,2,function(x,codon) as.matrix(codon)%*%as.matrix(x),codon)
# maxval = apply(CU,2,max)
# CU = apply(CU,1,function(x,maxval) x/maxval,maxval)
CU = t(apply(CU,2,AAnormalize,paste0(codons[rownames(codon),"AA"],rownames(codon)))) # normalize by AA

# Calculate normalized TE based on demand and supply of codons
nTE = array(NaN, c(ncol(anticodon), ncol(anticodon), nrow(codon)));
colnames(nTE) = colnames(anticodon); rownames(nTE) = colnames(anticodon); dimnames(nTE)[[3]] = rownames(codon)
initial_s = c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)

# Calculate tAI
for (sample in colnames(trna)){
  # Calculate relative adaptiveness values (ws)
  sample.ws = get.ws(tRNA=anticodon[,sample], s=initial_s, sking=0)
  sample.ws = AAnormalize(sample.ws,paste0(codons[rownames(codon),"AA"],rownames(codon))) # normalize by AA
  # sample.ws = sample.ws/max(sample.ws)
  
  # Calculate codon supply based on expression
  cTE = apply(CU,1,function(x,sample.ws) sample.ws/x,sample.ws)
  # maxval = apply(cTE,2,max)
  # nTE[,sample,] = apply(cTE,1,function(x,maxval) x/maxval,maxval)
  nTE[,sample,] = t(cTE) # keep non-normalized results
  # nTE[,sample,] = t(apply(cTE,2,AAnormalize,paste0(codons[rownames(codon),"AA"],rownames(codon)))) # normalize by AA
}

save(nTE,file="/home/xhernandez/Documents/tAI_genomic/AAnTE_CUgenomic_sqrt.rda")