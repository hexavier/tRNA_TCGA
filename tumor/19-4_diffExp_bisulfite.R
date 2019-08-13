library(ggplot2)
library(ggpubr)

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
anticodon = transformdata(trna,"")
# Upload results from differential analysis
diffexp = read.csv("results/differential_TvsH_twosided.csv", row.names = 1)
# Count delta catypes up vs down
delta = sapply(rownames(anticodon), function(x) (sum(diffexp[(diffexp$isoacceptor==x)&(diffexp$pval_adj<0.05),"FC"]>0,na.rm = T) - 
                                                   sum(diffexp[(diffexp$isoacceptor==x)&(diffexp$pval_adj<0.05),"FC"]<0,na.rm = T)))
bottom10 = names(delta)[delta<=quantile(delta,0.15625)]
top10 = names(delta)[delta>=quantile(delta,0.84375)]
notchanging = names(delta)[delta==0]

### METHYLATION ###
values = read.csv("results/differentialBisulfite_TvsH_twosided.csv",row.names = 1)
# Remove NA
values = values[!is.na(values$delta),]
# Classify genes into groups
values$trnagroup = NA
values[values$isoacceptor %in% top10,"trnagroup"] = "TOP"
values[values$isoacceptor %in% notchanging,"trnagroup"] = "NEUTRAL"
values[values$isoacceptor %in% bottom10,"trnagroup"] = "BOTTOM"
values = values[!is.na(values$trnagroup),]

# Plot
ggplot(values, aes(x=catype, y=delta, fill=trnagroup)) +  
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values=c("blue","white", "red")) +
  labs(title="Methylation across TCGA",x="tRNA group", y = "log2(FC) Methylation") + 
  stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
  stat_compare_means(aes(label = ..p.signif..)) +
  theme_classic()
