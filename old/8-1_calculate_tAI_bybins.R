## Load trna and weighted CU
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# tRNAs
trna = read.csv("data/TCGAall_nomod.csv",row.names = 1)
anticodon = extract_cod(transformdata(trna,"arcsinh"),codons$ANTICODON)

bins=20
for (n in 1:bins){
  # Codon usage
  weighted_CU = read.csv(sprintf("results/CUbybins/healthy_CU_bin%s.csv",(n-1)),row.names = 1)
  short_samples = sapply(colnames(trna), function(x) paste(strsplit(x,"\\.")[[1]][2:5],collapse="."))
  matched_CU = data.frame(sapply(short_samples, function(x) if(any((substr(colnames(weighted_CU),1,16) %in% x))){
    weighted_CU[,(substr(colnames(weighted_CU),1,16) %in% x)]
  }else{matrix(data=NA,nrow=1,ncol=64)}))
  rownames(matched_CU) = sapply(rownames(weighted_CU),function(x) paste(codons[x,"AA"],x,sep=""))
  codon = extract_cod(transformdata(matched_CU,""),rownames(codons)[!(codons$AA %in% c("Stop","Met"))])
  
  TAIs = data.frame(matrix(ncol = ncol(anticodon), nrow = ncol(anticodon)),row.names = colnames(anticodon)); colnames(TAIs) = colnames(anticodon)
  ## Calculate tAI
  for (sample in colnames(trna)){
    # Calculate relative adaptiveness values (ws)
    sample.ws = get.ws(tRNA=anticodon[,sample], sking=0)
    # Calculate tAI for all CUs
    sample.tai <- get.tai(t(codon), sample.ws)
    TAIs[sample,] = sample.tai
  }
  
  write.csv(TAIs,sprintf("results/CUbybins/tAI_relcodons_bin%s.csv",(n-1)))
}
