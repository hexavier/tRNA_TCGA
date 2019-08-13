library(ggplot2)
library(ggpubr)

mean_anticodon <- function(df){
  anticodons = sapply(rownames(df),function(x) strsplit(x,"-")[[1]][2])
  df_out = t(sapply(unique(anticodons), function(x) colMeans(df[anticodons==x,], na.rm=T)))
  return(df_out)
}

# Load data
cancer_types = c(BRCA="BRCA",PRAD="PRAD",kidney="KICH;KIRP;KIRC",lung="LUAD;LUSC",HNSC="HNSC",uterus="UCEC;CESC",
                 liver="LIHC;CHOL",THCA="THCA",colorectal="COAD;READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG")
path="/users/lserrano/xhernandez/tRNA_methylation/"

## Differential analysis Methylation
dataset = c()
for (t in names(cancer_types)){
  # Get data
  allvalues = read.csv(sprintf("%sResults/%s_trnas_TSS1500meth.csv",path,t), row.names = 1)
  types = unlist(strsplit(cancer_types[t],";"))
  for (type in types){
    # Get cancer type samples
    type_samples = colnames(read.csv(sprintf("%sData/Meth/%s.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt",path,type), sep="\t", nrows=1))
    type_samples = unique(substr(type_samples[2:length(type_samples)],1,15))
    # Subtract values
    rawvalues = allvalues[,(substr(colnames(allvalues),1,15) %in% type_samples)]
    # Average by anticodons
    meanvalues = mean_anticodon(rawvalues)
    
    ## Analyze methylation status of genes
    expr_samples = substr(colnames(meanvalues),1,15)
    mapped_normal = (regexpr("TCGA\\.[A-Z,0-9]{2}\\.[A-Z,0-9]{4}\\.11",expr_samples))==1
    
    dataset_temp = data.frame(row.names = colnames(meanvalues))
    dataset_temp$sample = expr_samples
    dataset_temp[,rownames(meanvalues)] = t(meanvalues)
    dataset_temp$catype = type
    dataset_temp$state = "CA"
    dataset_temp[mapped_normal,"state"] = "HE"
    dataset = rbind(dataset, dataset_temp)
  }
}

pdf("plots/Meth_diffexp_boxplots_rel.pdf",height=5, width=12)
for (gene in rownames(meanvalues)){
  print(ggplot(dataset, aes(x=catype, y=get(gene), fill=state)) +  
          geom_boxplot(position=position_dodge(1)) +
          scale_fill_manual(values=c("red", "blue")) +
          labs(title=sprintf("%s Methylation across TCGA",gene),x="Cancer Type", y = "Methylation") + 
          stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
          stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
          theme_classic())
}
dev.off()


## Differential analysis CNA
path="/users/lserrano/xhernandez/tRNA_scna/"
cancer_types = c(BRCA="BRCA",PRAD="PRAD",kidney="KICH;KIRP;KIRC",lung="LUAD;LUSC",HNSC="HNSC",uterus="UCEC;CESC",
                 liver="LIHC;CHOL",THCA="THCA",colorectal="COAD;READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG")
dataset = c()
for (t in names(cancer_types)){
  # Get data
  allvalues = read.csv(sprintf("%sResults/%s_trnas_cna.csv",path,t), row.names = 1)
  types = unlist(strsplit(cancer_types[t],";"))
  for (type in types){
    # Get cancer type samples
    type_samples = read.csv(sprintf("%sData/CNA/%s.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.seg.txt",path,type), sep="\t", colClasses = c("character",rep("NULL",5)))
    type_samples = gsub("-", ".", substr(unique(type_samples$Sample),1,15))
    # Subtract values
    rawvalues = allvalues[,(substr(colnames(allvalues),1,15) %in% type_samples)]
    # Average by anticodons
    meanvalues = mean_anticodon(rawvalues)
    
    ## Analyze methylation status of genes
    expr_samples = substr(colnames(meanvalues),1,15)
    mapped_normal = (regexpr("TCGA\\.[A-Z,0-9]{2}\\.[A-Z,0-9]{4}\\.11",expr_samples))==1
    
    dataset_temp = data.frame(row.names = colnames(meanvalues))
    dataset_temp$sample = expr_samples
    dataset_temp[,rownames(meanvalues)] = t(meanvalues)
    dataset_temp$catype = type
    dataset_temp$state = "CA"
    dataset_temp[mapped_normal,"state"] = "HE"
    dataset = rbind(dataset, dataset_temp)
  }
}

pdf("plots/CNA_diffexp_boxplots_rel.pdf",height=5, width=12)
for (gene in rownames(meanvalues)){
  print(ggplot(dataset, aes(x=catype, y=get(gene), fill=state)) +  
          geom_boxplot(position=position_dodge(1)) +
          scale_fill_manual(values=c("red", "blue")) +
          labs(title=sprintf("%s CNA across TCGA",gene),x="Cancer Type", y = "CNA") + 
          stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
          stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test") +
          theme_classic())
}
dev.off()