
#%% Load data
protnTE = read.csv("results/AAcTE_CUproteomics_sqrt.csv",row.names = 1)
seqnTE = read.csv("results/AAcTE_CUgenomic_sqrt.csv",row.names = 1); seqnTE = seqnTE[,!is.na(seqnTE[1,])]
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)
rownames(protnTE) = paste0(codons[rownames(protnTE),"AA"],rownames(protnTE))
rownames(seqnTE) = paste0(codons[rownames(seqnTE),"AA"],rownames(seqnTE))

# Correlate samples
samples = colnames(protnTE)[colnames(protnTE) %in% colnames(seqnTE)]
correlations = data.frame(sapply(samples,function(x) cor(protnTE[,x],seqnTE[,x],method = "spearman")))
colnames(correlations) = "Spearman Rho"
write.csv(correlations,"results/protTE_VS_seqTE_correlations.csv")

## Correlate expression
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
expcor = data.frame(SpearmanRho = numeric())
for (type in names(cancer_types)){
  # Get data
  seqdata = get_expression(cancer_types[type])
  colnames(seqdata) = substr(colnames(seqdata),1,16)
  protdata = read.csv(sprintf("/home/xhernandez/Downloads/TCGA-proteomics/TCGA_%s.csv",type),row.names = 1)
  # Find matching genes and samples
  genes = rownames(protdata)[rownames(protdata) %in% rownames(seqdata)]
  samples = colnames(protdata)[colnames(protdata) %in% colnames(seqdata)]
  # Correlate
  temp_cor = data.frame(sapply(samples,function(x) cor(protdata[genes,x],seqdata[genes,x],method = "spearman",use="na.or.complete")))
  rownames(temp_cor) = sprintf("%s-%s",type,samples)
  expcor = rbind(expcor,temp_cor)
}
colnames(expcor) = "Spearman Rho"
write.csv(expcor,"results/prot_VS_seq_correlations.csv")