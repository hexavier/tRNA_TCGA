library(eulerr)
library(gplots)

consistent_modif <- function(dataset, groups, thres=0.7){
  output = list()
  for (s in unique(groups)){
    idxs = groups %in% s
    tempdata = dataset[idxs,]
    counts = summary(tempdata$modif_id,maxsum=nrow(tempdata))
    # Keep only modifications in over 70% of samples
    nsamples = length(unique(tempdata$sample))
    output[[s]] = names(counts[counts>=(thres*nsamples)])
  }
  return(output)
}

mean_modif <- function(dataset, groups){
  allmodif = as.character(unique(dataset$modif_id))
  output = data.frame(row.names = allmodif)
  for (s in unique(groups)){
    idxs = groups %in% s
    tempdata = dataset[idxs,]
    output[allmodif,s] = sapply(allmodif, function(x) if(any(tempdata$modif_id %in% x)){mean(tempdata[(tempdata$modif_id %in% x),"percent"])}else{NA})
  }
  return(output)
}

## Load data
cancer_types = c(BRCA="BRCA",PRAD="PRAD",kidney="KICH;KIRP;KIRC",lung="LUAD;LUSC",HNSC="HNSC",uterus="UCEC;CESC",
                 liver="LIHC;CHOL",THCA="THCA",colorectal="COAD;READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG",GBM="GBM")
cancer_types = c(BRCA="BRCA",PRAD="PRAD",KICH="KICH",KIRC="KIRC",KIRP="KIRP",LUAD="LUAD",LUSC="LUSC",HNSC="HNSC",UCEC="UCEC",CESC="CESC",
                 LIHC="LIHC", CHOL="CHOL",THCA="THCA",COAD="COAD",READ="READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG",GBM="GBM")
rawmodH = read.csv("data/TCGAhealthy_modifications.csv",colClasses = c("NULL","character","integer",rep("character", 3),rep("integer", 8),"character","character"))
rawmodT = read.csv("data/TCGAtumor_modifications.csv",colClasses = c("NULL","character","integer",rep("character", 3),rep("integer", 8),"character","character"))
rawmod = rbind(rawmodH,rawmodT)

## Add trna info to cluster
clusters = as.character(read.csv("/users/lserrano/xhernandez/tRNA_mapping/Data/Genomes/H.sapiens/hg38.tRNAscan_clusterInfo.fa",
                                 header = 0)$V1)[c(T,F)]
mapclust = data.frame(cluster = as.character(sapply(clusters,function(x) strsplit(x,":|>")[[1]][2])),
                      trna = as.character(sapply(clusters,function(x) strsplit(x,"-|\\(")[[1]][2])))
mapclust = unique(mapclust); rownames(mapclust) = mapclust$cluster
rawmod$modif_id = as.factor(sapply(rawmod$modif_id, function(x) sprintf("%s-%s",mapclust[strsplit(as.character(x),"-")[[1]][1],"trna"],x)))

## Plot diagram
# HEK
fit <- euler(c(SMALL = sum(!(smallmod[["HEK293"]] %in% hydromod[["HEK"]])), 
               HYDRO = sum(!(hydromod[["HEK"]] %in% smallmod[["HEK293"]])), 
               "SMALL&HYDRO" = sum(smallmod[["HEK293"]] %in% hydromod[["HEK"]])))
print(plot(fit, fills=c("darkgoldenrod1","darkseagreen3"), edges=F, legend=F, labels=T, quantities=T, main="HEK"))

## Find consistent changes within replicates
tissue = sapply(as.character(rawmod$sample), function(x) strsplit(x,"-")[[1]][1])
tissue = sapply(tissue,function(y) names(cancer_types)[sapply(cancer_types, function(x) any(strsplit(x,";")[[1]] %in% y))])
tum_code = sapply(regexpr("[A-Z]{3,4}-TCGA-[A-Z,0-9]{2}-[A-Z,0-9]{4}-11",as.character(rawmod$sample)), function(x) if(x==1){"HE"}else{"CA"})
groups=sprintf("%s-%s",tissue,tum_code)
modif = consistent_modif(rawmod,groups,0.5)

##Save results
events = unique(as.character(unlist(modif)))
output = data.frame(row.names = events)
for (s in unique(groups)){
  output[,s] = sapply(events,function(x) x %in% modif[[s]])
}
write.csv(output,"results/modifications_tissues.csv")

# Plot heatmap
binarytab = apply(output,2,as.numeric); rownames(binarytab) = rownames(output)
binarytab = binarytab[rowSums(binarytab)>1,]
heatmap.2(binarytab,col=c("red", "blue"),symm=F, margins = c(6, 10))

# Plot VENN diagrams
pdf("plots/venn_modification_overlap.pdf")
## Plot diagram
for (type in names(cancer_types)){
  fit <- euler(c(TUMOR = sum(!(modif[[sprintf("%s-%s",type,"CA")]] %in% modif[[sprintf("%s-%s",type,"HE")]])), 
                 NORMAL = sum(!(modif[[sprintf("%s-%s",type,"HE")]] %in% modif[[sprintf("%s-%s",type,"CA")]])), 
                 "TUMOR&NORMAL" = sum(modif[[sprintf("%s-%s",type,"CA")]] %in% modif[[sprintf("%s-%s",type,"HE")]])))
  print(plot(fit, fills=c("darkgoldenrod1","darkseagreen3"), edges=F, legend=F, labels=T, quantities=T, main=type))
} 
dev.off()

###### Analyze A-to-I ##########
rawmod$percent = rawmod$altCount/rawmod$rawDepth
quant_modif = mean_modif(rawmod,groups)
modif34 = grep("[A-Za-z]{3}A[ACTG]{2}-cluster[0-9]+-3[0-9]-AtoG",rownames(quant_modif))
AtoImodif = quant_modif[modif34,]

# Plot heatmap
heatmap.2(as.matrix(AtoImodif), Rowv=F,symm=T, margins = c(4, 12))