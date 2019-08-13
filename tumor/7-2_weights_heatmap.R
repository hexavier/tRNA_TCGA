library(gplots)

# Load data
discdata = as.matrix(read.csv("results/AAcTE_discriminant_analyses.csv", row.names = 1))
# Codon table
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)
rownames(discdata) = paste0(codons[rownames(discdata),"AA"],rownames(discdata))

#Plot
breaks = c(min(discdata),seq(quantile(discdata,0.1),quantile(discdata,0.9),length.out=254),max(discdata))
heatmap.2(discdata[order((rowSums(discdata,na.rm=T)),decreasing = T),],
          Rowv=F, Colv=F, col = bluered(255), symm=T, margins=c(9,15), trace="none",
          ylab="delta PSI", xlab="Cancer Types", na.color="grey",breaks = breaks)
