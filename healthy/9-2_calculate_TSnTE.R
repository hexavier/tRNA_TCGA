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

#FROM ARTICLE: A benchmark of gene expression tissue-specificity metrics
#Function require a data frame with expression data, and give back a vector with PEM scores
#If expression for one tissue is not known, gene specificity for this gene is NA
fPem <- function(x)
{
  if(!all(is.na(x)))
  {
    x <- as.matrix(x)
    x[x<0] <- NA
    x <- cbind(x, r=rowSums(x, na.rm=FALSE)) #Add column with expression of gene per tissue
    x <- rbind(x, c=colSums(x, na.rm=TRUE))	#Add row with expression of all genes in a given tissue
    x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] <- x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] / (x[which(x[,ncol(x)]>0), ncol(x)] %o% x[nrow(x), which(x[nrow(x),]>0)] / x[nrow(x), ncol(x)]) #calculate the score
    
    x[x<1] <- 1
    x <- log10(x)
    
    x<- abs(x)				
    res <- x[-nrow(x),-ncol(x)] # remove the extra added columns
    res <- res/max(res, na.rm=TRUE) #Modification: To bring to normalized scale from 0 to 1
  } else {
    res <- NA
    print("No data avalable.")
  }
  return(res)
}

## Load trna and weighted CU
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# tRNAs
trna = read.csv("data/TCGAall_nomod.csv",row.names = 1)
# Rename columns
colnames(trna) = gsub("KICH","kidney",colnames(trna))
colnames(trna) = gsub("KIRP","kidney",colnames(trna))
colnames(trna) = gsub("KIRC","kidney",colnames(trna))
colnames(trna) = gsub("LUAD","lung",colnames(trna))
colnames(trna) = gsub("LUSC","lung",colnames(trna))
colnames(trna) = gsub("COAD","colorectal",colnames(trna))
colnames(trna) = gsub("READ","colorectal",colnames(trna))
colnames(trna) = gsub("CHOL","liver",colnames(trna))
colnames(trna) = gsub("LIHC","liver",colnames(trna))
colnames(trna) = gsub("UCEC","uterus",colnames(trna))
colnames(trna) = gsub("CESC","uterus",colnames(trna))

anticodon = extract_cod(transformdata(trna,"sqrt"),codons$ANTICODON)

# Genomic codon usage
codus = read.csv("data/human_CU_refseq.tsv",sep="\t")

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

## Calculate tissue-specificity PEM index
expr_bytissue = data.frame(row.names = rownames(isofseq))
tissues = sapply(colnames(anticodon), function(x) strsplit(x,"\\.")[[1]][1])
tissues_unique = unique(tissues)
# For each tissue, calculate Tau
for (t in tissues_unique){
  tissue_bool = (tissues==t)
  tissue_mean = rowMeans(matched_genexp[,tissue_bool],na.rm=T)
  expr_bytissue[,t] = tissue_mean
}
PEM = fPem(expr_bytissue)
write.csv(PEM,"results/TSclustering_index_pergene.csv")
PEM = read.csv("results/TSclustering_index_pergene.csv", row.names = 1)

###### SUPPLY CENTRIC ######
# Calculate normalized TE based on demand and supply of codons
nTE = array(NaN, c(ncol(anticodon), ncol(anticodon), nrow(codon)));
colnames(nTE) = colnames(anticodon); rownames(nTE) = colnames(anticodon); dimnames(nTE)[[3]] = rownames(codon)
initial_s = c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)

# Calculate tAI
for (sample in colnames(trna)){
  # Calculate for samples of that tissue
  t = strsplit(sample,"\\.")[[1]][1]
  # Identify and select tissue specific genes
  topTS = PEM[,t]>=0.5
  # Calculate CU of tissue specific genes
  CU = apply(matched_genexp[topTS,],2,function(x,codon) as.matrix(codon)%*%as.matrix(x),codon[,topTS])
  CU = t(apply(CU,2,AAnormalize,paste0(codons[rownames(codon),"AA"],rownames(codon)))) # normalize by AA
  
  # Calculate relative adaptiveness values (ws)
  sample.ws = get.ws(tRNA=anticodon[,sample], s=initial_s, sking=0)
  sample.ws = AAnormalize(sample.ws,paste0(codons[rownames(codon),"AA"],rownames(codon))) # normalize by AA

  # Calculate codon supply based on expression
  cTE = apply(CU,1,function(x,sample.ws) sample.ws/x,sample.ws)
  nTE[,sample,] = t(cTE) # keep non-normalized results
}

save(nTE,file="/home/xhernandez/Documents/tAI_genomic/AAcTE_TSCUgenomic_sqrt_SUPcentric.rda")

###### DEMAND CENTRIC ######

# Calculate normalized TE based on demand and supply of codons
nTE = array(NaN, c(ncol(anticodon), ncol(anticodon), nrow(codon)));
colnames(nTE) = colnames(anticodon); rownames(nTE) = colnames(anticodon); dimnames(nTE)[[3]] = rownames(codon)
initial_s = c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)

## Calculate tAI for tissue specific CU
tissues = sapply(colnames(trna), function(x) strsplit(x,"\\.")[[1]][1])
tissues_unique = unique(tissues)
# Calculate tAI
for (sample in colnames(trna)){
  for (t in tissues_unique){
    # Calculate for samples of that tissue
    idx_tissue = (tissues==t)
    # Identify and select tissue specific genes
    topTS = PEM[,t]>=0.5
    
    # Calculate CU of tissue specific genes
    CU = apply(matched_genexp[topTS,idx_tissue],2,function(x,codon) as.matrix(codon)%*%as.matrix(x),codon[,topTS])
    CU = t(apply(CU,2,AAnormalize,paste0(codons[rownames(codon),"AA"],rownames(codon)))) # normalize by AA
    
    # Calculate relative adaptiveness values (ws)
    sample.ws = get.ws(tRNA=anticodon[,sample], s=initial_s, sking=0)
    sample.ws = AAnormalize(sample.ws,paste0(codons[rownames(codon),"AA"],rownames(codon))) # normalize by AA
    
    # Calculate codon supply based on expression
    cTE = apply(CU,1,function(x,sample.ws) sample.ws/x,sample.ws)
    nTE[idx_tissue,sample,] = t(cTE) # keep non-normalized results
  }
}

save(nTE,file="/home/xhernandez/Documents/tAI_genomic/AAcTE_TSCUgenomic_sqrt_DEMcentric.rda")