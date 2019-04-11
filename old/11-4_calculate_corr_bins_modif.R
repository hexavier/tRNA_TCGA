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
TAIs = read.csv("/home/xhernandez/Downloads/tAI_genomic/tAI_modif_CUgenomic.csv", row.names = 1)
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

# Binarize data
bins = 50
items_per_bin = nrow(TAIs)/bins
# Initiate data structure
corTAI_CU = data.frame(matrix(ncol = ncol(TAIs), nrow = ncol(TAIs)),row.names = colnames(TAIs)); colnames(corTAI_CU) = colnames(TAIs)
corTAI_trna = data.frame(matrix(ncol = ncol(TAIs), nrow = ncol(TAIs)),row.names = colnames(TAIs)); colnames(corTAI_trna) = colnames(TAIs)

for (s in colnames(TAIs)){
  temp_expr = sapply(1:bins, function(x) colMeans(asinh(matched_genexp[order(TAIs[,s])[trunc(items_per_bin*(x-1)):trunc(items_per_bin*x)],]), na.rm = T))
  temp_tai = sapply(1:bins, function(x) colMeans(TAIs[order(TAIs[,s])[trunc(items_per_bin*(x-1)):trunc(items_per_bin*x)],], na.rm = T))

  # Calculate CU
  corTAI_trna[,s] = cor(t(temp_expr),temp_tai[s,],method="spearman")
  corTAI_CU[s,] = t(cor(temp_expr[s,],t(temp_tai),method="spearman"))
}

# Calculate weighted tAI
write.csv(corTAI_trna,"results/binscor_modif_tAI_trnacentric.csv")
write.csv(corTAI_CU,"results/binscor_modif_tAI_CUcentric.csv")
