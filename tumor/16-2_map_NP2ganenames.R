# load data
TAIs = read.csv("results/deltaRTE_KIRC.csv",row.names = 1)

# Load mapping file
np2nm = read.csv("data/NP2NM.txt",sep="\t")
#http://genome-euro.ucsc.edu/cgi-bin/hgTables ncbiRefLink
np2nm = np2nm[np2nm[,6]!="",]
rownames(np2nm) = sapply(np2nm[,6], function(x) substr(as.character(x),1,nchar(as.character(x))-2))

# Rename
proteins = as.character(sapply(rownames(TAIs), function(x) substr(as.character(x),1,nchar(as.character(x))-2)))
genes = as.character(np2nm[proteins,"name"])

# Remove unmapped NPs
TAIs = TAIs[!is.na(genes),]
genes=genes[!is.na(genes)]

# Find duplicate genes and average them
ugenes = unique(genes)
TAIs_out = t(sapply(ugenes,function(x) colMeans(TAIs[genes==x,],na.rm = T)))

# Save
write.csv(TAIs_out,"results/deltaRTE_KIRC_renamed.csv")