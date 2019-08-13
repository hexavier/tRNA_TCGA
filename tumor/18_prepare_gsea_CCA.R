transformdata <- function(data,transf){
  aa_idx = regexpr("i?[A-Z][a-z]{2}[A-Z]{3}",rownames(data))==1
  data = data[aa_idx,]
  if (transf=="log"){
    outdata = sapply(data,log)
    # Remove inf values
    outdata[outdata==-Inf] = NaN
    rownames(outdata)=rownames(data)
  }else if (transf=="arcsinh"){
    outdata = sapply(data,asinh)
    rownames(outdata)=rownames(data)
  }else if (transf=="sqrt"){
    outdata = sapply(data,sqrt)
    rownames(outdata)=rownames(data)
  }else if (transf=="rel"){
    # Compute relative data
    outdata = data.frame(matrix(ncol = ncol(data), nrow = nrow(data)),row.names = rownames(data))
    colnames(outdata)= colnames(data)
    aa = sapply(rownames(outdata),function(x) substr(x,1,nchar(x)-3))
    uniqueaa = unique(aa)
    for (n in uniqueaa){
      idx = (aa %in% n)
      idx_data = matrix(as.matrix(data[idx,]), ncol = ncol(data), nrow = sum(idx))
      total = colSums(idx_data)
      outdata[idx,] = t(apply(idx_data,1,function(x) x/total))
      iszero = (total %in% 0)
      if (any(iszero)){
        outdata[idx,iszero] = 1.0/sum(idx)
      }
    }
  }else{
    outdata=data
  }
  return(outdata)
}

# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# Genomic codon usage
codus = read.csv("data/human_CU_refseq.tsv",sep="\t")
codus_idx = as.character(codus$Protein.ID)
# Keep only columns with codon info
codus_clean = data.frame(sapply(unique(codus_idx),function(x) colMeans(codus[(codus_idx==x),13:ncol(codus)],na.rm = T)))
rownames(codus_clean) = sapply(rownames(codus_clean),function(x) paste(codons[x,"AA"],x,sep=""))

# Retrieve codons
codon = transformdata(codus_clean,"rel")
# Keep CCA
codon = t(codon["ProCCA",])

## Rename proteins
# Load mapping file
np2nm = read.csv("data/NP2NM.txt",sep="\t")
#http://genome-euro.ucsc.edu/cgi-bin/hgTables ncbiRefLink
np2nm = np2nm[np2nm[,6]!="",]
rownames(np2nm) = sapply(np2nm[,6], function(x) substr(as.character(x),1,nchar(as.character(x))-2))

# Rename
proteins = as.character(sapply(rownames(codon), function(x) substr(as.character(x),1,nchar(as.character(x))-2)))
genes = as.character(np2nm[proteins,"name"])

# Remove unmapped NPs
codon = codon[!is.na(genes),]
genes=genes[!is.na(genes)]

# Find duplicate genes and average them
ugenes = unique(genes)
codon_out = t(sapply(ugenes,function(x) mean(codon[genes %in% x],na.rm = T)))

# Save
write.csv(t(codon_out),"results/relCCAcontent_genes.csv")