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
    mapped_tumor = (regexpr("TCGA\\.[A-Z,0-9]{2}\\.[A-Z,0-9]{4}\\.01",expr_samples))==1
    
    dataset_temp = data.frame(row.names = colnames(meanvalues))
    dataset_temp$sample = expr_samples
    dataset_temp[,rownames(meanvalues)] = t(meanvalues)
    dataset_temp$catype = type
    dataset_temp$state = NA
    dataset_temp[mapped_normal,"state"] = "HE"
    dataset_temp[mapped_tumor,"state"] = "CA"
    dataset_temp = dataset_temp[!is.na(dataset_temp$state),]
    dataset = rbind(dataset, dataset_temp)
  }
}

# Build summary table statistics
cancer_types = c(BRCA="BRCA",PRAD="PRAD",KICH="KICH",KIRC="KIRC",KIRP="KIRP",LUAD="LUAD",LUSC="LUSC",HNSC="HNSC",UCEC="UCEC",CESC="CESC",
                 LIHC="LIHC", CHOL="CHOL",THCA="THCA",COAD="COAD",READ="READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG")

summary_table = data.frame(row.names = 1:(length(cancer_types)*nrow(meanvalues)))
catypes = c()
gene = c()
meansT = c()
meansH = c()
pvals = c()
for (g in rownames(meanvalues)){
  catypes = append(catypes,names(cancer_types))
  gene = append(gene, rep(g,length(cancer_types)))
  meansT = append(meansT,sapply(names(cancer_types),
                                function(x,dataset) mean(dataset[(dataset$catype==x)&(dataset$state=="CA"),g],na.rm=T),dataset))
  meansH = append(meansH,sapply(names(cancer_types),
                                function(x,dataset) mean(dataset[(dataset$catype==x)&(dataset$state=="HE"),g],na.rm=T),dataset))
  # Calculate pvalue if at least 2 healthy and 2 tumor samples are available
  pvals = append(pvals,sapply(names(cancer_types),
                              function(x,dataset) if (all(c(sum(!is.na(dataset[(dataset$catype==x)&(dataset$state=="HE"),g])),
                                                            sum(!is.na(dataset[(dataset$catype==x)&(dataset$state=="CA"),g])))>1)){
                                wilcox.test(dataset[(dataset$catype==x)&(dataset$state=="HE"),g],
                                            dataset[(dataset$catype==x)&(dataset$state=="CA"),g],
                                            alternative = "two.sided")$p.value}else{NA},dataset))
}
summary_table$catype = catypes
summary_table$gene = gene
summary_table$mean_exp_T = meansT
summary_table$mean_exp_H = meansH
summary_table$FC = log2(meansT/meansH)
summary_table$pval = pvals
summary_table$pval_adj = p.adjust(pvals,method="fdr")
write.csv(summary_table, sprintf("results/differentialMeanMeth_TvsH_twosided.csv",type))


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
    mapped_tumor = (regexpr("TCGA\\.[A-Z,0-9]{2}\\.[A-Z,0-9]{4}\\.01",expr_samples))==1
    
    dataset_temp = data.frame(row.names = colnames(meanvalues))
    dataset_temp$sample = expr_samples
    dataset_temp[,rownames(meanvalues)] = t(meanvalues)
    dataset_temp$catype = type
    dataset_temp$state = NA
    dataset_temp[mapped_normal,"state"] = "HE"
    dataset_temp[mapped_tumor,"state"] = "CA"
    dataset_temp = dataset_temp[!is.na(dataset_temp$state),]
    dataset = rbind(dataset, dataset_temp)
  }
}

# Build summary table statistics
cancer_types = c(BRCA="BRCA",PRAD="PRAD",KICH="KICH",KIRC="KIRC",KIRP="KIRP",LUAD="LUAD",LUSC="LUSC",HNSC="HNSC",UCEC="UCEC",CESC="CESC",
                 LIHC="LIHC", CHOL="CHOL",THCA="THCA",COAD="COAD",READ="READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG")

summary_table = data.frame(row.names = 1:(length(cancer_types)*nrow(meanvalues)))
catypes = c()
gene=c()
meansT = c()
meansH = c()
pvals = c()
for (g in rownames(meanvalues)){
  catypes = append(catypes,names(cancer_types))
  gene = append(gene, rep(g,length(cancer_types)))
  meansT = append(meansT,sapply(names(cancer_types),
                                function(x,dataset) mean(dataset[(dataset$catype==x)&(dataset$state=="CA"),g],na.rm=T),dataset))
  meansH = append(meansH,sapply(names(cancer_types),
                                function(x,dataset) mean(dataset[(dataset$catype==x)&(dataset$state=="HE"),g],na.rm=T),dataset))
  # Calculate pvalue if at least 2 healthy and 2 tumor samples are available
  pvals = append(pvals,sapply(names(cancer_types),
                              function(x,dataset) if (all(c(sum(!is.na(dataset[(dataset$catype==x)&(dataset$state=="HE"),g])),
                                                            sum(!is.na(dataset[(dataset$catype==x)&(dataset$state=="CA"),g])))>1)){
                                wilcox.test(dataset[(dataset$catype==x)&(dataset$state=="HE"),g],
                                            dataset[(dataset$catype==x)&(dataset$state=="CA"),g],
                                            alternative = "two.sided")$p.value}else{NA},dataset))
}
summary_table$catype = catypes
summary_table$gene = gene
summary_table$mean_exp_T = meansT
summary_table$mean_exp_H = meansH
summary_table$delta = (meansT - meansH)
summary_table$pval = pvals
summary_table$pval_adj = p.adjust(pvals,method="fdr")
write.csv(summary_table, sprintf("results/differentialMeanCNA_TvsH_twosided.csv",type))
