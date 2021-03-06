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
TAIs = read.csv("/home/xhernandez/Documents/tAI_genomic/tAI_CUgenomic_sqrt.csv", row.names = 1)
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

# Calculate weighted tAI
weightedTAI = as.matrix(t(matched_genexp))%*%as.matrix(TAIs)
weightedTAI = weightedTAI/colSums(matched_genexp)
write.csv(weightedTAI,"results/weighted_tAI.csv")

# Calculate correlations
corTAI = data.frame(matrix(ncol = ncol(TAIs), nrow = ncol(TAIs)),row.names = colnames(TAIs)); colnames(corTAI) = colnames(TAIs)

for (s in colnames(TAIs)){
  corTAI[,s] = cor(matched_genexp,TAIs[,s],method="spearman")
}
write.csv(corTAI,"results/cor_tAI.csv")

