library(gplots)
plotdots <- function(s,n){
  nrows = round(sqrt(n))
  # Define point vectors
  xdots = unlist(sapply(1:nrows, function(x) if(n/x>=nrows){rep(x,nrows)}else{rep(x,(n%%nrows))}))
  ydots = c(rep(1:nrows,trunc(sqrt(n))),nrows:(nrows-(n%%nrows)+1))
  # Separate successes and failures
  x_suc = xdots[1:s]
  y_suc = ydots[1:s]
  x_fail = xdots[(s+1):length(xdots)]
  y_fail = ydots[(s+1):length(ydots)]
  # Plot scatter
  plot(xdots,ydots,type="n",axes=FALSE, frame.plot=FALSE,ann=FALSE)
  points(x_suc,y_suc,pch = 19, col="blue",cex=4)
  points(x_fail,y_fail,pch = 19, col="red",cex=4)
}

# Load MLR coefficients
meth = read.csv("results/MLR_methylation_coefficients.csv", row.names = 1)
cna = read.csv("results/MLR_cna_coefficients.csv", row.names = 1)

# Analyze whether cna affects significantly positively
x = sum(cna>0, na.rm = T)
n = sum(!is.na(cna))
binom.test(x,n)
# Plot points
plotdots(x,n)

# Analyze whether cna affects significantly positively
x = sum(meth<0, na.rm = T)
n = sum(!is.na(meth))
binom.test(x,n)
# Plot points
plotdots(n-x,n)

### PLOT HEATMAP
diffexp_mean = as.matrix(meth)
# Heatmap
breaks = c(min(diffexp_mean,na.rm=T),seq(-50,50,length.out=254),max(diffexp_mean,na.rm=T))
diffexp_mean = diffexp_mean[order((rowSums(diffexp_mean<0,na.rm=T)-rowSums(diffexp_mean>0,na.rm=T)),decreasing = T),order(colSums(is.na(diffexp_mean)))]
heatmap.2(diffexp_mean,
          Rowv=F, Colv=F, col = bluered(255), symm=T, margins=c(6,6), trace="none",
          na.color="grey",breaks = breaks)

### PLOT Bloxplot
diffexp_mean = as.matrix(meth)
diffexp_mean = diffexp_mean[order((rowSums(diffexp_mean<0,na.rm=T)-rowSums(diffexp_mean>0,na.rm=T)),decreasing = T),order(colSums(is.na(diffexp_mean)))]

# Re-structure data
dataset = data.frame(row.names = 1:(ncol(diffexp_mean)*nrow(diffexp_mean)))
dataset$value = as.numeric(diffexp_mean)
dataset$trna = as.character(sapply(rownames(diffexp_mean), function(x) rep(x,ncol(diffexp_mean))))
dataset$catype = as.character(rep(colnames(diffexp_mean), nrow(diffexp_mean)))
dataset = dataset[!is.na(dataset$value),]
# Plot
ggplot(dataset, aes(x=trna, y=value)) +  
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values=c("black")) +
  labs(title="Coefficients across TCGA",x="tRNA", y = "Coefficient") + 
  stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) + 
  theme_classic()

