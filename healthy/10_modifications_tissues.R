library(UpSetR)
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
rawmod = read.csv("data/TCGAall_modifications.csv")

## Add trna info to cluster
clusters = as.character(read.csv("/users/lserrano/xhernandez/tRNA_mapping/Data/Genomes/H.sapiens/hg38.tRNAscan_clusterInfo.fa",
                                 header = 0)$V1)[c(T,F)]
mapclust = data.frame(cluster = as.character(sapply(clusters,function(x) strsplit(x,":|>")[[1]][2])),
                      trna = as.character(sapply(clusters,function(x) strsplit(x,"-|\\(")[[1]][2])))
mapclust = unique(mapclust); rownames(mapclust) = mapclust$cluster
rawmod$modif_id = as.factor(sapply(rawmod$modif_id, function(x) sprintf("%s-%s",mapclust[strsplit(as.character(x),"-")[[1]][1],"trna"],x)))

## Find consistent changes within replicates
groups = sapply(as.character(rawmod$sample), function(x) strsplit(x,"-")[[1]][1])
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
# Keep only modification detected in at least 2 tissues
binarytab = binarytab[rowSums(binarytab)>1,]
heatmap.2(binarytab,col=c("red", "blue"),symm=F, margins = c(4, 10))

###### Analyze A-to-I ##########
rawmod$percent = rawmod$altCount/rawmod$rawDepth
quant_modif = mean_modif(rawmod,groups)
modif34 = grep("[A-Za-z]{3}A[ACTG]{2}-cluster[0-9]+-3[0-9]-AtoG",rownames(quant_modif))
AtoImodif = quant_modif[modif34,]

# Plot heatmap
heatmap.2(as.matrix(AtoImodif), Rowv=F,symm=T, margins = c(4, 12))