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
trnaH = read.csv("data/TCGAhealthy_nomod.csv",row.names = 1)
trnaT = read.csv("data/TCGAtumor_nomod.csv",row.names = 1)
trna = cbind(trnaH,trnaT)
anticodon = transformdata(trna,"sqrt")

# Average by anticodons
rawvalues = read.csv(sprintf("%sResults/trnas_bisulfite.csv",path), row.names = 1)
#rawvalues = read.csv(sprintf("%sResults/trnas_TSS1500bisulfite.csv",path), row.names = 1)
values = mean_anticodon(rawvalues)
# Match samples
trna_samples = substr(sapply(colnames(anticodon), function(x) paste(strsplit(x,"\\.")[[1]][4:5],collapse=".")),1,7)
value_samples = sapply(colnames(values), function(x) if(nchar(strsplit(x,"_")[[1]][3])==4){paste(strsplit(x,"_")[[1]][3],"01",sep = ".")}
                       else{paste(substr(strsplit(x,"_")[[1]][3],2,5),"11",sep = ".")})
merged_trna = colnames(anticodon)[trna_samples %in% value_samples]
merged_value_idx = sapply(trna_samples[trna_samples %in% value_samples], function(x) which(x==value_samples))
# Match anticodons
acods = rownames(values)[rownames(values) %in% rownames(anticodon)]
# Calculate correlation with expression
correlations = data.frame(sapply(acods,function(x) cor(anticodon[x,merged_trna],
                                                   values[x,merged_value_idx],method = "spearman",use="na.or.complete")))
colnames(correlations) = "Spearman Methylation"

# Plot
hist(as.numeric(correlations$`Spearman Methylation`))