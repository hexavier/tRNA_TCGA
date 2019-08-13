library(RTCGAToolbox)

## Define cancer types
cancer_abbr = c("BRCA","PRAD","KIRC","KICH","KIRP","LUAD","LUSC","HNSC","UCEC","LIHC","THCA",
                "COAD","READ","PAAD","STAD","CHOL","SKCM","THYM","ESCA","BLCA","GBM","CESC","PCPG")

## Calculate mitotic index
for (type in cancer_abbr){
  dwload = getFirehoseData(dataset=type, destdir="/home/xhernandez/Downloads/TCGA-mRNAseq", RNASeq2GeneNorm=TRUE)
  data = getData(dwload,"RNASeq2GeneNorm")
  subset = (regexpr("TCGA-[A-Z,0-9]{2}-[A-Z,0-9]{4}-11",colnames(data)))==1
  if (!exists("healthy")){
    healthy = data[,subset]
  }else{
    healthy = cbind(healthy,data[,subset])
  }
}

write.csv(healthy,file="data/healthy_mrnaseq.csv")