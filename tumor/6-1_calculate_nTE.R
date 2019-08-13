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

get_expression <- function(abbr){
  if (sum(grep(";",abbr))>0){
    names = unlist(strsplit(abbr,";"))
    for (n in names){
      if (!exists("output")){
        output = read.csv(sprintf("/home/xhernandez/Downloads/TCGA-isoformSeq/%s.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.data.txt",n), 
                          row.names = 1, sep= "\t",na.strings="normalized_count")
        output = output[2:nrow(output),]
      }else{
        toadd = read.csv(sprintf("/home/xhernandez/Downloads/TCGA-isoformSeq/%s.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.data.txt",n), 
                         row.names = 1, sep= "\t",na.strings="normalized_count")
        toadd = toadd[2:nrow(toadd),]
        output = cbind(output, toadd)
      }
    }
  }else{
    output = read.csv(sprintf("/home/xhernandez/Downloads/TCGA-isoformSeq/%s.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.data.txt",abbr), 
                      row.names = 1, sep= "\t",na.strings="normalized_count")
    output = output[2:nrow(output),]
  }
  return(output)
}

# Define cancer types
cancer_types = c(BRCA="BRCA",PRAD="PRAD",kidney="KICH;KIRP;KIRC",lung="LUAD;LUSC",HNSC="HNSC",uterus="UCEC;CESC",
                 liver="LIHC;CHOL",THCA="THCA",colorectal="COAD;READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG")

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
tcga2refseq = read.csv("data/TCGA2RefSeq.txt",sep="\t", header=F,row.names = 1)
# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownIsoforms.txt.gz
rownames(tcga2refseq) = sapply(rownames(tcga2refseq), function(x) substr(x,1,nchar(x)-2))

np2nm = read.csv("data/NP2NM.txt",sep="\t", header=T)
#http://genome-euro.ucsc.edu/cgi-bin/hgTables ncbiRefLink
np2nm = np2nm[np2nm[,6]!="",]
rownames(np2nm) = sapply(np2nm[,6], function(x) substr(as.character(x),1,nchar(as.character(x))-2))

# Isoform gene expression
isofseq = read.csv("data/healthy_isoform_seq.csv",row.names = 1)
rownames(isofseq) = as.character(sapply(rownames(isofseq), function(x) substr(x,1,nchar(x)-2)))

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

# Prepare codon data
rownames(codus_clean) = sapply(rownames(codus_clean),function(x) paste(codons[x,"AA"],x,sep=""))
codon = extract_cod(transformdata(codus_clean,""),rownames(codons)[!(codons$AA %in% c("Stop","Met"))])

# Calculate normalized TE based on demand and supply of codons
nTE = matrix(NaN,ncol = ncol(anticodon), nrow = nrow(codon));
colnames(nTE) = colnames(anticodon); rownames(nTE) = rownames(codon)
initial_s = c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)

for (type in names(cancer_types)){
  idxtissue = (sapply(colnames(anticodon), function(x) strsplit(x,"\\.")[[1]][1]) %in% strsplit(cancer_types[type],";")[[1]])
  # Substract gene expression and trna data
  temptrna = anticodon[,idxtissue]

  # Clean isofseq dataset
  isofseq = get_expression(cancer_types[type])
  rownames(isofseq) = as.character(sapply(rownames(isofseq), function(x) substr(x,1,nchar(x)-2)))
  isofseq = isofseq[tcga_ids,]
  
  # Map gene expression samples
  short_samples = sapply(colnames(temptrna), function(x) paste(strsplit(x,"\\.")[[1]][2:5],collapse="."))
  short_genexp_samples = substr(colnames(isofseq),1,16)
  matched_genexp = data.frame(sapply(short_samples, function(x) if(any((short_genexp_samples==x))){
    isofseq[,(short_genexp_samples==x)]}else{matrix(data=NA,nrow=1,ncol=ncol(codon))}))
  CU = apply(matched_genexp,2,function(x,codon) as.matrix(codon)%*%as.matrix(x),codon)
  CU = t(apply(CU,2,AAnormalize,paste0(codons[rownames(codon),"AA"],rownames(codon)))) # normalize by AA
  
  # Calculate tAI
  sample.ws = apply(temptrna,2, get.ws, s=initial_s, sking=0)
  sample.ws = t(apply(sample.ws,2,AAnormalize,paste0(codons[rownames(codon),"AA"],rownames(codon)))) # normalize by AA
  
  # Compute ratio between supply and demand
  cTE = sample.ws/CU
  nTE[,idxtissue] = t(cTE) # keep non-normalized results
}

write.csv(nTE,"results/AAcTE_CUgenomic_sqrt.csv")