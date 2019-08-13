# tRNAs
trna = read.csv("data/TCGAall_nomod.csv",row.names = 1)

# Genomic codon usage
codus = read.csv("data/human_CU_refseq.tsv",sep="\t")

# Isoform gene expression
isofseq = read.csv("data/healthy_isoform_seq.csv",row.names = 1)
rownames(isofseq) = as.character(sapply(rownames(isofseq), function(x) substr(x,1,nchar(x)-2)))

# Mapping
tcga2refseq = read.csv("data/TCGA2RefSeq.txt",sep="\t", header=F,row.names = 1)
# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownIsoforms.txt.gz
rownames(tcga2refseq) = sapply(rownames(tcga2refseq), function(x) substr(x,1,nchar(x)-2))

## Prepare CU data
# Map ids
tcga_ids = rownames(isofseq)[(rownames(isofseq) %in% rownames(tcga2refseq))]
refseq_ids = sapply(rownames(isofseq), function(x) if (x %in% rownames(tcga2refseq)){substr(tcga2refseq[x,1],3,100)}else{NA})
refseq_ids = as.character(refseq_ids[!is.na(refseq_ids)])
# Detect whether refseq ids are in codon usage data, and remove if not
codusindex_main = sapply(codus$Protein.ID, function(x) substr(x,3,nchar(as.character(x))-2))
# Keep only events in common between tcga and codon usage
keep = (refseq_ids %in% codusindex_main)
tcga_ids = tcga_ids[keep]
refseq_ids = refseq_ids[keep]
# Convert booleans in codus_ids to string index. If more than 1, take mean
codus_clean = data.frame(sapply(refseq_ids,function(x) colMeans(codus[(codusindex_main==x),13:ncol(codus)],na.rm = T)))

# Clean isofseq dataset
isofseq = isofseq[tcga_ids,]

# Match samples
short_samples = sapply(colnames(trna), function(x) paste(strsplit(x,"\\.")[[1]][2:5],collapse="."))
short_genexp_samples = substr(colnames(isofseq),1,16)
matched_genexp = data.frame(sapply(short_samples, function(x) if(any((short_genexp_samples==x))){
  isofseq[,(short_genexp_samples==x)]}else{matrix(data=NA,nrow=1,ncol=nrow(isofseq))}))

## Clustering
expr=matched_genexp
expr_bytissue = data.frame(row.names = rownames(isofseq))
tissues = sapply(colnames(expr), function(x) strsplit(x,"\\.")[[1]][1])
tissues_unique = unique(tissues)
# For each tissue, calculate expression
for (t in tissues_unique){
  tissue_bool = (tissues==t)
  tissue_mean = rowMeans(expr[,tissue_bool],na.rm=T)
  expr_bytissue[,t] = tissue_mean
}

clusters <- hclust(dist(t(expr_bytissue),method="manhattan"))
plot(clusters)
clusterCut <- cutree(clusters, 17) # establish cut at 17 clusters
table(clusterCut, tissues_unique)

