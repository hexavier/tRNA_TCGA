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


## Calculate tAI for genomic CU
rownames(codus_clean) = sapply(rownames(codus_clean),function(x) paste(codons[x,"AA"],x,sep=""))
codon = extract_cod(transformdata(codus_clean,""),rownames(codons)[!(codons$AA %in% c("Stop","Met"))])

TAIs = data.frame(matrix(ncol = ncol(anticodon), nrow = ncol(codon)),row.names = colnames(codon)); colnames(TAIs) = colnames(anticodon)
initial_s = c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)
# Calculate tAI
for (sample in colnames(trna)){
  # Calculate relative adaptiveness values (ws)
  sample.ws = get.ws(tRNA=anticodon[,sample], s=initial_s, sking=0)
  #sample.ws = AAnormalize(sample.ws,paste0(codons[rownames(codon),"AA"],rownames(codon))) # normalize by AA
  # Calculate tAI for all CUs
  sample.tai <- get.tai(t(codon), sample.ws)
  TAIs[,sample] = sample.tai
}

write.csv(TAIs,"/home/xhernandez/Documents/tAI_genomic/tAI_CUgenomic_sqrt.csv")

TAIs = read.csv("/home/xhernandez/Documents/tAI_genomic/tAI_CUgenomic_sqrt.csv", row.names = 1)

## Analyze correlation with gene expression
# Map samples
short_samples = sapply(colnames(TAIs), function(x) paste(strsplit(x,"\\.")[[1]][2:5],collapse="."))
short_genexp_samples = substr(colnames(isofseq),1,16)
matched_genexp = data.frame(sapply(short_samples, function(x) if(any((short_genexp_samples==x))){
  isofseq[,(short_genexp_samples==x)]}else{matrix(data=NA,nrow=1,ncol=nrow(TAIs))}))

# Remove 0s
#nozero = rowSums(matched_genexp==0,na.rm = T)==0
#matched_genexp = matched_genexp[nozero,]
#TAIs = TAIs[nozero,]

# Calculate whether tAI correlates with expression
#cor_bysample = diag(cor(TAIs,asinh(matched_genexp),method="spearman"))

bins = 50
items_per_bin = nrow(TAIs)/bins
cors=c()
for (s in colnames(TAIs)){
  #bins_expr = sapply(1:bins, function(x) mean(sort(asinh(matched_genexp[,s]))[trunc(items_per_bin*(x-1)):trunc(items_per_bin*x)]))
  #bins_tai = sapply(1:bins, function(x) mean(TAIs[order(matched_genexp[,s]),s][trunc(items_per_bin*(x-1)):trunc(items_per_bin*x)]))
  #plot(bins_expr,bins_tai)
  bins_expr = sapply(1:bins, function(x) mean(asinh(matched_genexp[order(TAIs[,s]),s])[trunc(items_per_bin*(x-1)):trunc(items_per_bin*x)]))
  bins_tai = sapply(1:bins, function(x) mean(sort(TAIs[,s])[trunc(items_per_bin*(x-1)):trunc(items_per_bin*x)]))
  ## Check how standart deviation change and plot error bars
  # sd_expr = sapply(1:bins, function(x) sd(asinh(matched_genexp[order(TAIs[,s]),s])[trunc(items_per_bin*(x-1)):trunc(items_per_bin*x)]))
  # sd_tai = sapply(1:bins, function(x) sd(sort(TAIs[,s])[trunc(items_per_bin*(x-1)):trunc(items_per_bin*x)]))
  # plot(bins_expr,bins_tai)
  # arrows(x0=bins_expr, y0=(bins_tai-sd_tai), x1=bins_expr, y1=(bins_tai+sd_tai),
  #        angle=90, code=3, length=0.04, lwd=0.4)
  # arrows(x0=(bins_expr-sd_expr), y0=bins_tai, x1=(bins_expr+sd_expr), y1=bins_tai,
  #        angle=90, code=3, length=0.04, lwd=0.4)
  ## Check a linear regression of all points
  # linreg = lm(tai ~ exp,data=data.frame(exp=asinh(matched_genexp[,s]),tai=TAIs[,s]))
  # print(linreg)
  c = cor(bins_expr,bins_tai,method="spearman")
  cors = append(cors,c)
}

tissues = sapply(colnames(TAIs), function(x) strsplit(x,"\\.")[[1]][1])
tissues_unique = unique(tissues)
#cor_bytissue = sapply(tissues_unique, function(x) median(cor_bysample[tissues %in% x],na.rm = T))
corbin_bytissue = sapply(tissues_unique, function(x) median(cors[tissues %in% x],na.rm = T))
dataset = data.frame(cors,tissues)
ggplot(dataset, aes(x=tissues, y=cors, fill=tissues)) +  
  geom_violin(position=position_dodge(1)) +
  labs(title="CU 50 bins correlation by tissue",x="Tissue", y = "Spearman Rho tAI vs Expr") + 
  stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
  theme_classic()
