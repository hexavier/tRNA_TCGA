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
TAIs = read.csv("results/tAI_rel_CUgenomic.csv", row.names = 1)

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

## Calculate tAI for top expressed genes
weightedTAI_cu = data.frame(matrix(ncol = ncol(TAIs), nrow = ncol(TAIs)),row.names = colnames(TAIs)); colnames(weightedTAI_cu) = colnames(TAIs)
weightedTAI_trna = data.frame(matrix(ncol = ncol(TAIs), nrow = ncol(TAIs)),row.names = colnames(TAIs)); colnames(weightedTAI_trna) = colnames(TAIs)

# Calculate tAI TOP
for (sample in colnames(TAIs)){
  for (s in colnames(TAIs)){
    # identify top expressed genes
    top = order(matched_genexp[,s], decreasing = T)[1:1000]
    # Identify and select tissue specific genes
    exprTS = matched_genexp[top,s]
    taiTS = TAIs[top,sample]
    # Calculate CU
    weightedTAI_cu[s,sample] = (as.matrix(t(exprTS))%*%as.matrix(taiTS))/sum(exprTS)
    
    # Identify and select tissue specific genes
    exprTS = matched_genexp[top,sample]
    taiTS = TAIs[top,s]
    # Calculate CU
    weightedTAI_trna[sample,s] = (as.matrix(t(exprTS))%*%as.matrix(taiTS))/sum(exprTS)
  }
}

write.csv(weightedTAI_cu,"results/TOPweighted_tAI_CUcentric.csv")
write.csv(weightedTAI_trna,"results/TOPweighted_tAI_trnacentric.csv")

## Calculate tAI for bottom expressed genes
weightedTAI_cu = data.frame(matrix(ncol = ncol(TAIs), nrow = ncol(TAIs)),row.names = colnames(TAIs)); colnames(weightedTAI_cu) = colnames(TAIs)
weightedTAI_trna = data.frame(matrix(ncol = ncol(TAIs), nrow = ncol(TAIs)),row.names = colnames(TAIs)); colnames(weightedTAI_trna) = colnames(TAIs)

# Calculate tAI BOTTOM
for (sample in colnames(TAIs)){
  for (s in colnames(TAIs)){
    # identify top expressed genes
    exprTS = matched_genexp[matched_genexp[,s]!=0,s]
    taiTS = TAIs[matched_genexp[,s]!=0,sample]
    top = order(exprTS, decreasing = F)[1:1000]
    # Identify and select tissue specific genes
    exprTS = exprTS[top]
    taiTS = taiTS[top]
    # Calculate CU
    weightedTAI_cu[s,sample] = (as.matrix(t(exprTS))%*%as.matrix(taiTS))/sum(exprTS)
    
    # identify top expressed genes
    exprTS = matched_genexp[matched_genexp[,sample]!=0,sample]
    taiTS = TAIs[matched_genexp[,sample]!=0,s]
    top = order(exprTS, decreasing = F)[1:1000]
    # Identify and select tissue specific genes
    exprTS = exprTS[top]
    taiTS = taiTS[top]
    # Calculate CU
    weightedTAI_trna[sample,s] = (as.matrix(t(exprTS))%*%as.matrix(taiTS))/sum(exprTS)
  }
}

write.csv(weightedTAI_cu,"results/BOTTOMweighted_tAI_CUcentric.csv")
write.csv(weightedTAI_trna,"results/BOTTOMweighted_tAI_trnacentric.csv")
