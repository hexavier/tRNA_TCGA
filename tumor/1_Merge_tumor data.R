# Define cancer types
cancer_types = c("BRCA","PRAD","kidney","lung","HNSC","uterus","liver","THCA","colorectal","ESCA","STAD","BLCA","PAAD","THYM",
  "SKCM","PCPG")
# Merge all quantification data
for (catype in cancer_types){
  if (exists("output")){
    temp = read.csv(sprintf("TCGA_%s/postprocessing/results_nomod.csv", catype),row.names = 1)
    output = cbind(output,temp)
  }else{
    output = read.csv(sprintf("TCGA_%s/postprocessing/results_nomod.csv", catype),row.names = 1)
  }
}
colnames(output)= gsub("\\.","-",colnames(output))
# Write data
write.csv(output, "/home/xhernandez/Dropbox/PhD/tRNA_TCGAtumor/data/TCGAtumor_nomod.csv")

# Merge all modification data
for (catype in cancer_types){
  if (exists("modif")){
    temp = read.csv(sprintf("TCGA_%s/postprocessing/modifications.csv", catype),row.names = 1)
    modif = rbind(modif,temp)
  }else{
    modif = read.csv(sprintf("TCGA_%s/postprocessing/modifications.csv", catype),row.names = 1)
  }
}

# Write data
write.csv(modif, "/home/xhernandez/Dropbox/PhD/tRNA_TCGAtumor/data/TCGAtumor_modifications.csv")