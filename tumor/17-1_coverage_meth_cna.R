change_chr <- function(chromosome){
  chr = strsplit(chromosome,"_|\\.|-")[[1]][1]
  chrnumb = substr(chr,4,5)
  return(chrnumb)
}

# Load data
cancer_types = c(BRCA="BRCA",PRAD="PRAD",kidney="KICH;KIRP;KIRC",lung="LUAD;LUSC",HNSC="HNSC",uterus="UCEC;CESC",
                 liver="LIHC;CHOL",THCA="THCA",colorectal="COAD;READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG")
# tRNA genes
path="/users/lserrano/xhernandez/tRNA_methylation/"
trna_coord = read.csv(paste0(path,"Data/Genomes/H.sapiens/hg19.tRNAscan.bed12"), sep="\t", header = F, row.names = 4)
trna_coord$V1 = sapply(as.character(trna_coord$V1),change_chr)

# Keep unique coordinates
genes = rownames(unique(trna_coord[,c(1,2,3)]))

## How many genes detected?
percentmeth = c()
percentcna = c()
for (t in names(cancer_types)){
  # Get data
  path="/users/lserrano/xhernandez/tRNA_methylation/"
  methraw = read.csv(sprintf("%sResults/%s_trnas_TSS1500meth.csv",path,t), row.names = 1)
  path="/users/lserrano/xhernandez/tRNA_scna/"
  cnaraw = read.csv(sprintf("%sResults/%s_trnas_cna.csv",path,t), row.names = 1)
  # Calculate %
  addmeth = apply(methraw,2,function(x) sum(!is.na(x))/length(x))
  percentmeth =  append(percentmeth, addmeth)
  
  addcna = apply(cnaraw,2,function(x) sum(!is.na(x))/length(x))
  percentcna =  append(percentcna, addcna)
}