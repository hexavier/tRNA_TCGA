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
# tRNA genes
path="/users/lserrano/xhernandez/tRNA_methylation/"
trna_coord = read.csv(paste0(path,"Data/Genomes/H.sapiens/hg19.tRNAscan.bed12"), sep="\t", header = F, row.names = 4)
trna_coord$V1 = sapply(as.character(trna_coord$V1),change_chr)
# tRNAs
trnaH = read.csv("data/TCGAhealthy_nomod.csv",row.names = 1)
trnaT = read.csv("data/TCGAtumor_nomod.csv",row.names = 1)
trna = cbind(trnaH,trnaT)
anticodon = transformdata(trna,"sqrt")

# Keep unique coordinates
genes = rownames(unique(trna_coord[,c(1,2,3)]))

## ANCOVA
ancovas= list()
for (t in names(cancer_types)){
  # Get data
  path="/users/lserrano/xhernandez/tRNA_methylation/"
  methvalues = read.csv(sprintf("%sResults/%s_trnas_TSS1500meth.csv",path,t), row.names = 1)
  path="/users/lserrano/xhernandez/tRNA_scna/"
  cnavalues = read.csv(sprintf("%sResults/%s_trnas_cna.csv",path,t), row.names = 1)
  types = unlist(strsplit(cancer_types[t],";"))
  for (type in types){
    # Get cancer type samples
    path="/users/lserrano/xhernandez/tRNA_methylation/"
    meth_samples = colnames(read.csv(sprintf("%sData/Meth/%s.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt",path,type), sep="\t", nrows=1))
    meth_samples = unique(substr(meth_samples[2:length(meth_samples)],1,16))
    path="/users/lserrano/xhernandez/tRNA_scna/"
    cna_samples = read.csv(sprintf("%sData/CNA/%s.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.seg.txt",path,type), sep="\t", colClasses = c("character",rep("NULL",5)))
    cna_samples = gsub("-", ".", substr(unique(cna_samples$Sample),1,16))
    # Subtract values
    value_samples = meth_samples[meth_samples %in% cna_samples]
    methraw = methvalues[,(substr(colnames(methvalues),1,16) %in% value_samples)]
    cnaraw = cnavalues[,(substr(colnames(cnavalues),1,16) %in% value_samples)]
    
    # Average by anticodons
    methmean = mean_anticodon(methraw)
    cnamean = mean_anticodon(cnaraw)
    
    # Match samples
    trna_samples = sapply(colnames(anticodon), function(x) paste(strsplit(x,"\\.")[[1]][2:5],collapse="."))
    meth_samples = substr(colnames(methmean),1,16)
    cna_samples = substr(colnames(cnamean),1,16)
    merged_trna = colnames(anticodon)[trna_samples %in% value_samples]
    merged_meth_idx = sapply(trna_samples[trna_samples %in% meth_samples], function(x) which(x==meth_samples))
    merged_cna_idx = sapply(trna_samples[trna_samples %in% cna_samples], function(x) which(x==cna_samples))
    
    # Match anticodons
    acods = rownames(methmean)[rownames(methmean) %in% rownames(anticodon)]
    # Build dataset
    dataset_temp = data.frame(row.names = 1:(length(merged_trna)*length(acods)))
    dataset_temp$sample = as.character(sapply(merged_trna,function(x) rep(x,length(acods))))
    dataset_temp$codon = as.factor(rep(acods, length(merged_trna)))
    dataset_temp$catype = as.factor(type)
    dataset_temp$methylation = as.numeric(sapply(merged_meth_idx,function(x) methmean[acods,x]))
    dataset_temp$cna = as.numeric(sapply(merged_cna_idx,function(x) cnamean[acods,x]))
    dataset_temp$expression = as.numeric(sapply(merged_trna,function(x) anticodon[acods,x]))
    # Remove rows that are NA or expression==0
    dataset_temp = dataset_temp[!is.na(dataset_temp[,"expression"]),]
    
    # Run ANCOVA
    ancovas[[type]] = lm(data = dataset_temp, expression ~ methylation + cna + codon + codon:methylation + codon:cna)
  }
}