
## Define cancer types
cancer_abbr = c("BRCA","PRAD","KIRC","KICH","KIRP","LUAD","LUSC","HNSC","UCEC","LIHC","THCA",
                "COAD","READ","PAAD","STAD","CHOL","SKCM","THYM","ESCA","BLCA","GBM","CESC","PCPG")

## Calculate mitotic index
for (type in cancer_abbr){
  path="/home/xhernandez/Downloads/TCGA-isoformSeq"
  filename = sprintf("%s.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.data.txt",type)
  data = read.delim(file.path(path,filename),header = T, sep = "\t", dec = ".",row.names = 1)
  colnam = gsub("\\.","-",colnames(data))
  subset = (regexpr("TCGA-[A-Z,0-9]{2}-[A-Z,0-9]{4}-11",colnam))==1
  data = data.frame(data[,subset]) 
  colnames(data) = colnam[subset]
  if (!exists("healthy")){
    healthy = data
  }else{
    healthy = cbind(healthy,data)
  }
}

write.csv(healthy,file="data/healthy_isoform_seq.csv")