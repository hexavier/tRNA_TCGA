change_chr <- function(chromosome){
  chr = strsplit(chromosome,"_|\\.|-")[[1]][1]
  chrnumb = substr(chr,4,5)
  return(chrnumb)
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
mean_anticodon <- function(df){
  anticodons = sapply(rownames(df),function(x) strsplit(x,"-")[[1]][2])
  df_out = t(sapply(unique(anticodons), function(x) colMeans(df[anticodons==x,], na.rm=T)))
  return(df_out)
}

# Load data
cancer_types = c(BRCA="BRCA",PRAD="PRAD",kidney="KICH;KIRP;KIRC",lung="LUAD;LUSC",HNSC="HNSC",uterus="UCEC;CESC",
                 liver="LIHC;CHOL",THCA="THCA",colorectal="COAD;READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG")
path="/users/lserrano/xhernandez/tRNA_methylation/"
trna_coord = read.csv(paste0(path,"Data/Genomes/H.sapiens/hg19.tRNAscan.bed12"), sep="\t", header = F, row.names = 4)
trna_coord$V1 = sapply(as.character(trna_coord$V1),change_chr)
# tRNAs
trna = read.csv("data/TCGAall_nomod.csv",row.names = 1)
anticodon = transformdata(trna,"sqrt")

# Keep unique coordinates
genes = rownames(unique(trna_coord[,c(1,2,3)]))

## Analyze methylation
rm(rawvalues)
for (type in names(cancer_types)){
  # Get data
  if (exists("rawvalues")){
    add = read.csv(sprintf("%sResults/%s_trnas_TSS1500meth.csv",path,type), row.names = 1)
    rawvalues = cbind(rawvalues,add)
  }else{
    rawvalues = read.csv(sprintf("%sResults/%s_trnas_TSS1500meth.csv",path,type), row.names = 1)
  }
}

# Average by anticodons
rawvalues = rawvalues[genes,]
values = mean_anticodon(rawvalues)
# Match samples
trna_samples = sapply(colnames(anticodon), function(x) paste(strsplit(x,"\\.")[[1]][2:5],collapse="."))
value_samples = substr(colnames(values),1,16)
merged_trna = colnames(anticodon)[trna_samples %in% value_samples]
merged_value_idx = sapply(trna_samples[trna_samples %in% value_samples], function(x) which(x==value_samples))
# Match anticodons
acods = rownames(values)[rownames(values) %in% rownames(anticodon)]
# Calculate correlation with expression
correlations = data.frame(sapply(acods,function(x) cor(anticodon[x,merged_trna],
                                                       values[x,merged_value_idx],method = "spearman",use="na.or.complete")))
colnames(correlations) = "Spearman Methylation"
write.csv(correlations,"results/trnaH_methylation_corr_allsamples.csv")

## Analyze CNA
path="/users/lserrano/xhernandez/tRNA_scna/"
rm(rawvalues)
for (type in names(cancer_types)){
  # Get data
  if (exists("rawvalues")){
    add = read.csv(sprintf("%sResults/%s_trnas_cna.csv",path,type), row.names = 1)
    rawvalues = cbind(rawvalues,add)
  }else{
    rawvalues = read.csv(sprintf("%sResults/%s_trnas_cna.csv",path,type), row.names = 1)
  }
}

# Average by anticodons
rawvalues = rawvalues[genes,]
values = mean_anticodon(rawvalues)
# Match samples
trna_samples = sapply(colnames(anticodon), function(x) paste(strsplit(x,"\\.")[[1]][2:5],collapse="."))
value_samples = substr(colnames(values),1,16)
merged_trna = colnames(anticodon)[trna_samples %in% value_samples]
merged_value_idx = sapply(trna_samples[trna_samples %in% value_samples], function(x) which(x==value_samples))
# Match anticodons
acods = rownames(values)[rownames(values) %in% rownames(anticodon)]
# Calculate correlation with expression
correlations = data.frame(sapply(acods,function(x) cor(anticodon[x,merged_trna],
                                                       values[x,merged_value_idx],method = "spearman",use="na.or.complete")))

colnames(correlations) = "Spearman CNA"
write.csv(correlations,"results/trnaH_CNA_corr_allsamples.csv")