library(RTCGAToolbox)

get_expression <- function(abbr)
{
  if (sum(grep(";",abbr))>0){
    names = unlist(strsplit(abbr,";"))
    for (n in names){
      dwload = getFirehoseData(dataset=n, RNASeq2GeneNorm=TRUE, destdir="/home/xhernandez/Downloads/TCGA-mRNAseq")
      if (!exists("output")){
        output = getData(dwload,"RNASeq2GeneNorm")
      }else{
        toadd = getData(dwload,"RNASeq2GeneNorm")
        output = cbind(output, toadd)
      }
    }
  }else{
    dwload = getFirehoseData(dataset=abbr, RNASeq2GeneNorm=TRUE, destdir="/home/xhernandez/Downloads/TCGA-mRNAseq")
    output = getData(dwload,"RNASeq2GeneNorm")
  }
  return(output)
}

# tRNAs
tissues = c("COAD",	"READ","GBM","BRCA","STAD","KIRP","KICH","PRAD","HNSC","THCA","KIRC","BLCA",
            "LUSC","LUAD","LIHC","PAAD","UCEC","PCPG","SKCM","CHOL","ESCA","THYM","CESC")
genexp = get_expression(tissues[1])

# Compute mean of tissue
genexp_mean = data.frame(row.names = rownames(genexp))
for (t in unique(tissues)){
  genexp = get_expression(t)
  state = substr(sapply(colnames(genexp), function(x) strsplit(x, "-")[[1]][4]),1,2)
  genexp_mean[,sprintf("%s-HE",t)] = rowMeans(genexp[,(state=="11"),drop=F],na.rm = T)
  genexp_mean[,sprintf("%s-CA",t)] = rowMeans(genexp[,(state=="01"),drop=F],na.rm = T)
}

write.csv(genexp_mean,"results/GeneExp_catypes.csv")