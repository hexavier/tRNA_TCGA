
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

## Load TAIs and exprression
# Mapping
tcga2refseq = read.csv("data/TCGA2RefSeq.txt",sep="\t", header=F,row.names = 1)
# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownIsoforms.txt.gz
rownames(tcga2refseq) = sapply(rownames(tcga2refseq), function(x) substr(x,1,nchar(x)-2))

np2nm = read.csv("data/NP2NM.txt",sep="\t", header=T)
# http://genome-euro.ucsc.edu/cgi-bin/hgTables ncbiRefLink
np2nm = np2nm[np2nm[,6]!="",]
rownames(np2nm) = sapply(np2nm[,6], function(x) substr(as.character(x),1,nchar(as.character(x))-2))

# Genomic codon usage
codus = read.csv("data/human_CU_refseq.tsv",sep="\t")

# Isoform gene expression
isofseq = read.csv("data/healthy_isoform_seq.csv",row.names = 1)
rownames(isofseq) = as.character(sapply(rownames(isofseq), function(x) substr(x,1,nchar(x)-2)))
tcga_ids = rownames(isofseq)[(rownames(isofseq) %in% rownames(tcga2refseq))]
refseq_ids = sapply(rownames(isofseq), function(x) if (x %in% rownames(tcga2refseq)){as.character(tcga2refseq[x,1])}else{NA})
refseq_ids = as.character(refseq_ids[!is.na(refseq_ids)])
# Detect whether refseq ids are in codon usage data, and remove if not
codusindex_main = sapply(codus$Protein.ID, function(x) substr(as.character(np2nm[substr(x,1,nchar(as.character(x))-2),"X.id"]),1,nchar(as.character(x))-2))

keep = (refseq_ids %in% codusindex_main)
tcga_ids = tcga_ids[keep]
isofseq = isofseq[tcga_ids,]

# TAIs
TAIs = read.csv("/home/xhernandez/Downloads/tAI_genomic/tAI_CUgenomic.csv", row.names = 1)

# Rename samples based on clustering
colnames(TAIs) = gsub("KICH","kidney",colnames(TAIs))
colnames(TAIs) = gsub("KIRP","kidney",colnames(TAIs))
colnames(TAIs) = gsub("KIRC","kidney",colnames(TAIs))
colnames(TAIs) = gsub("LUAD","lung",colnames(TAIs))
colnames(TAIs) = gsub("LUSC","lung",colnames(TAIs))
colnames(TAIs) = gsub("COAD","colorectal",colnames(TAIs))
colnames(TAIs) = gsub("READ","colorectal",colnames(TAIs))
colnames(TAIs) = gsub("CHOL","liver",colnames(TAIs))
colnames(TAIs) = gsub("LIHC","liver",colnames(TAIs))
colnames(TAIs) = gsub("UCEC","uterus",colnames(TAIs))
colnames(TAIs) = gsub("CESC","uterus",colnames(TAIs))

# Map samples
short_samples = sapply(colnames(TAIs), function(x) paste(strsplit(x,"\\.")[[1]][2:5],collapse="."))
short_genexp_samples = substr(colnames(isofseq),1,16)
matched_genexp = data.frame(sapply(short_samples, function(x) if(any((short_genexp_samples==x))){
  isofseq[,(short_genexp_samples==x)]}else{matrix(data=NA,nrow=1,ncol=nrow(TAIs))}))

## Calculate tissue-specificity PEM index
expr_bytissue = data.frame(row.names = rownames(isofseq))
tissues = sapply(colnames(TAIs), function(x) strsplit(x,"\\.")[[1]][1])
tissues_unique = unique(tissues)
# For each tissue, calculate Tau
for (t in tissues_unique){
  tissue_bool = (tissues==t)
  tissue_mean = rowMeans(matched_genexp[,tissue_bool],na.rm=T)
  expr_bytissue[,t] = tissue_mean
}
PEM = fPem(expr_bytissue)
# write.csv(PEM,"results/TSclustering_index_pergene.csv")
# PEM = read.csv("results/TSclustering_index_pergene.csv", row.names = 1)


## Calculate tAI for tissue specific CU
# Initiate data structure
corTAI_CU = data.frame(matrix(ncol = ncol(TAIs), nrow = ncol(TAIs)),row.names = colnames(TAIs)); colnames(corTAI_CU) = colnames(TAIs)
corTAI_trna = data.frame(matrix(ncol = ncol(TAIs), nrow = ncol(TAIs)),row.names = colnames(TAIs)); colnames(corTAI_trna) = colnames(TAIs)

# Calculate tAI for CU centric
for (sample in colnames(TAIs)){
  # Calculate for samples of that tissue
  t = strsplit(sample,"\\.")[[1]][1]
  # Identify and select tissue specific genes
  topTS = PEM[,t]>=0.5
  exprTS = matched_genexp[topTS,sample]
  taiTS = TAIs[topTS,]

  # Calculate CU
  corTAI_CU[sample,] = t(cor(exprTS,taiTS,method="spearman"))
}

# Calculate tAI for trna centric
for (s in colnames(TAIs)){
  # Calculate for samples of that tissue
  t = strsplit(s,"\\.")[[1]][1]
  # Identify and select tissue specific genes
  topTS = PEM[,t]>=0.5
  exprTS = matched_genexp[topTS,]
  taiTS = TAIs[topTS,s]

  # Calculate CU
  corTAI_trna[,s] = cor(exprTS,taiTS,method="spearman")
}

write.csv(corTAI_CU,"results/TScor_tAI_CUcentric.csv")
write.csv(corTAI_trna,"results/TScor_tAI_trnacentric.csv")

## Calculate tAI for MOST VARIABLY EXPRESSED GENES
# Calculate coefficient of variation
coeffvar = apply(expr_bytissue,1, function(x) sd(x)/mean(x))
coeffvar[is.na(coeffvar)]=0

# Initiate data structure
corTAI = data.frame(matrix(ncol = ncol(TAIs), nrow = ncol(TAIs)),row.names = colnames(TAIs)); colnames(corTAI) = colnames(TAIs)

# Calculate tAI correlation
for (s in colnames(TAIs)){
  corTAI[,s] = cor(matched_genexp[coeffvar>2,],TAIs[coeffvar>2,s],method="spearman")
}
write.csv(corTAI,"results/corMOSTVAR_tAI.csv")