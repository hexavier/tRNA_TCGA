## Correlate expression with tAI
get_expression <- function(abbr){
  if (sum(grep(";",abbr))>0){
    names = unlist(strsplit(abbr,";"))
    for (n in names){
      if (!exists("output")){
        output = read.csv(sprintf("/home/xhernandez/Downloads/TCGA-mRNAseq/20160128-%s-RNAseq2GeneNorm.txt",n), 
                          row.names = 1, sep= "\t",na.strings="normalized_count")
        output = output[2:nrow(output),]
      }else{
        toadd = read.csv(sprintf("/home/xhernandez/Downloads/TCGA-mRNAseq/20160128-%s-RNAseq2GeneNorm.txt",n), 
                         row.names = 1, sep= "\t",na.strings="normalized_count")
        toadd = toadd[2:nrow(toadd),]
        output = cbind(output, toadd)
      }
    }
  }else{
    output = read.csv(sprintf("/home/xhernandez/Downloads/TCGA-mRNAseq/20160128-%s-RNAseq2GeneNorm.txt",abbr), 
                      row.names = 1, sep= "\t",na.strings="normalized_count")
    output = output[2:nrow(output),]
  }
  rnames=as.character(sapply(rownames(output),function(x) strsplit(x,"\\|")[[1]][1]))
  output = output[!duplicated(rnames),]; rnames=rnames[!duplicated(rnames)]
  rownames(output) = rnames
  return(output)
}

cancer_types = c(BRCA="BRCA",colorectal="COAD;READ")
correlations = data.frame(row.names = c("RTE_ratio","RTE_prot","tAI_ratio","tAI_prot","RtAI_ratio","RtAI_prot"))
for (type in names(cancer_types)){
  # Get data
  seqdata = get_expression(cancer_types[type])
  colnames(seqdata) = substr(colnames(seqdata),1,16)
  protdata = read.csv(sprintf("/home/xhernandez/Downloads/TCGA-proteomics/TCGA_%s.csv",type),row.names = 1)
  # Load tAIs
  tAI = read.csv(sprintf("results/tAI_CUproteomics%s_sqrt.csv",type),row.names = 1)
  colnames(tAI)=sapply(colnames(tAI), function(x) paste(strsplit(x,"\\.")[[1]][2:5],collapse = "."))
  RTE = read.csv(sprintf("results/RTEAI_CUproteomics%s_sqrt.csv",type),row.names = 1)
  colnames(RTE)=sapply(colnames(RTE), function(x) paste(strsplit(x,"\\.")[[1]][2:5],collapse = "."))
  RtAI = read.csv(sprintf("results/AAtAI_CUproteomics%s_sqrt.csv",type),row.names = 1)
  colnames(RtAI)=sapply(colnames(RtAI), function(x) paste(strsplit(x,"\\.")[[1]][2:5],collapse = "."))
  # Find matching genes and samples
  genes = rownames(tAI)[rownames(tAI) %in% rownames(protdata)[rownames(protdata) %in% rownames(seqdata)]]
  samples = colnames(tAI)[colnames(tAI) %in% colnames(protdata)[colnames(protdata) %in% colnames(seqdata)]]
  # Calculate ratios
  ratio = data.frame(sapply(samples,function(x) protdata[genes,x]/seqdata[genes,x]))
  prots = protdata[genes,samples]
  # Subset tAI and RTE
  RTE = RTE[genes,(colnames(RTE) %in% samples)]
  tAI = tAI[genes,(colnames(tAI) %in% samples)]
  RtAI = RtAI[genes,(colnames(RtAI) %in% samples)]
  # Calculate correlations
  correlations[,sprintf("%s-%s",type,samples)] = sapply(1:length(samples),function(x) as.numeric(cor(cbind(ratio[,x],prots[,x],RTE[,x],tAI[,x],RtAI[,x]),method="spearman",use="na.or.complete")[1:2,3:5]))
}
correlations=t(correlations)
write.csv(correlations,"results/tAIvsProt_correlations_bysample.csv")

## Plot correlations
library(ggplot2)
library(ggpubr)
# Structure data
dataset = c()
for (s in colnames(correlations)){
  dataset_temp = data.frame(row.names=seq(1,nrow(correlations)))
  dataset_temp$corr = as.numeric(correlations[,s])
  dataset_temp$sample = rownames(correlations)
  dataset_temp$TEmetric = strsplit(s,"_")[[1]][1]
  dataset_temp$expression = strsplit(s,"_")[[1]][2]
  dataset = rbind(dataset, dataset_temp)
}

# Plot
compare_means(corr ~ TEmetric,data=dataset,ref.group = "RTE", group.by = "expression",method="wilcox.test",paired = T)
compare_means(corr ~ expression,data=dataset,group.by = "TEmetric",method="wilcox.test",paired = T)

pdf("plots/tAI_proteomics_corr_bysample_paired.pdf",height=5,width = 9)
print(ggplot(dataset, aes(x=expression, y=corr, fill=TEmetric)) +  
        geom_violin(position=position_dodge(1)) +
        labs(title="TE Metric Comparison",x="Expression", y = "Spearman Correlation") + 
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
        stat_compare_means(aes(label = ..p.format..),method="wilcox.test",paired = T) +
        theme_classic())
dev.off()