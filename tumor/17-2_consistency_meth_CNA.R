library(gplots)

#%% Load data
cancer_types = c(BRCA="BRCA",PRAD="PRAD",KICH="KICH",KIRC="KIRC",KIRP="KIRP",LUAD="LUAD",LUSC="LUSC",HNSC="HNSC",UCEC="UCEC",CESC="CESC",
                 LIHC="LIHC", CHOL="CHOL",THCA="THCA",COAD="COAD",READ="READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG")
values = read.csv("results/differentialMeth_TvsH_twosided.csv",row.names = 1)
#values = read.csv("results/differentialCNA_TvsH_twosided.csv",row.names = 1)

#Initialize structure
anticodons = unique(values$isoacceptor)
consistency = matrix(nrow = length(anticodons), ncol=length(cancer_types))
rownames(consistency) = anticodons; colnames(consistency) = cancer_types
for (c in anticodons){
  aadata = values[values$isoacceptor==c,]
  up = sapply(cancer_types,function(x) sum(aadata[aadata$catype==x,"FC"]>0,na.rm = T)/sum(!is.na(aadata[aadata$catype==x,"FC"])))
  down = sapply(cancer_types,function(x) sum(aadata[aadata$catype==x,"FC"]<0,na.rm = T)/sum(!is.na(aadata[aadata$catype==x,"FC"])))
  consistency[c,] = (up-down)
}

# Heatmap
breaks = c(min(consistency,na.rm=T),seq(-0.5,0.5,length.out=254),max(consistency,na.rm=T))
diffexp_mean = consistency[order((rowSums(consistency<0,na.rm=T)-rowSums(consistency>0,na.rm=T)),decreasing = T),order(colSums(is.na(consistency)))]
heatmap.2(diffexp_mean,
          Rowv=F, Colv=F, col = bluered(255), symm=T, margins=c(12,15), trace="none",
          na.color="grey",breaks = breaks)
