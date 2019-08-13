library(survminer)
library(survival)
library(RTCGAToolbox)

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

binarize <- function(expr,patients_expr,patients_cli,qntl=0.5)
{
  mapped = sapply(patients_cli,function(x) grep(sprintf("%s.01",x),patients_expr))
  # Take average of patients with multiple samples
  cont_tum = as.numeric(lapply(mapped,function(x) mean(expr[x],na.rm=T)))
  bin = rep(NA, length(patients_cli))
  bin[cont_tum>quantile(cont_tum,1-qntl,na.rm=T)] = T
  bin[cont_tum<=quantile(cont_tum,qntl,na.rm=T)] = F
  return(bin)
}

get_clinical <- function(abbr)
{
  if (sum(grep(";",abbr))>0){
    names = unlist(strsplit(abbr,";"))
    for (n in names){
      dwload = getFirehoseData(dataset=n, clinical=TRUE, destdir = "/home/xhernandez/Downloads/TCGA-clinical")
      if (!exists("output")){
        output = getData(dwload,"clinical")
      }else{
        toadd = getData(dwload,"clinical")
        output = rbind(output[,c("vital_status","days_to_death","days_to_last_followup")],toadd[,c("vital_status","days_to_death","days_to_last_followup")])
      }
    }
  }else{
    dwload = getFirehoseData(dataset=abbr, clinical=TRUE, destdir = "/home/xhernandez/Downloads/TCGA-clinical")
    output = getData(dwload,"clinical")
  }
  output[!is.na(output$days_to_death),"times"] = output[!is.na(output$days_to_death),"days_to_death"]
  output[!is.na(output$days_to_last_followup),"times"] = output[!is.na(output$days_to_last_followup),"days_to_last_followup"]
  cleanoutput = data.frame(times=as.numeric(output$times),vital_status=as.numeric(output$vital_status),row.names = rownames(output))
  return(cleanoutput)
}


## Upload data

cancer_types = c(BRCA="BRCA",PRAD="PRAD",KICH="KICH",KIRC="KIRC",KIRP="KIRP",LUAD="LUAD",LUSC="LUSC",HNSC="HNSC",UCEC="UCEC",CESC="CESC",
                 LIHC="LIHC", CHOL="CHOL",THCA="THCA",COAD="COAD",READ="READ",ESCA="ESCA",STAD="STAD",BLCA="BLCA",PAAD="PAAD",THYM="THYM",
                 SKCM="SKCM",PCPG="PCPG")
trnaH = read.csv("data/TCGAhealthy_nomod.csv",row.names = 1)
trnaT = read.csv("data/TCGAtumor_nomod.csv",row.names = 1)
trna = cbind(trnaH,trnaT)
anticodon = transformdata(trna,"rel")
# Upload results from differential analysis
diffexp = read.csv("results/differentialREL_TvsH_twosided.csv", row.names = 1)
# Count delta catypes up vs down
delta = sapply(rownames(anticodon), function(x) (sum(diffexp[(diffexp$isoacceptor==x)&(diffexp$pval_adj<0.05),"FC"]>0,na.rm = T) - 
                                                   sum(diffexp[(diffexp$isoacceptor==x)&(diffexp$pval_adj<0.05),"FC"]<0,na.rm = T)))
filtered_trnas = names(delta[abs(delta)>=6])

# Create summary table
summary = data.frame(row.names = names(cancer_types))
for (type in names(cancer_types)){
  # Clinical data
  survInfo = get_clinical(cancer_types[type])
  
  # Retrieve RTE data
  idxtissue = (sapply(colnames(anticodon), function(x) strsplit(x,"\\.")[[1]][1]) %in% strsplit(cancer_types[type],";")[[1]])
  tempdata = anticodon[,idxtissue]
  colnames(tempdata) = sapply(colnames(tempdata), function(x) paste(strsplit(x,"\\.")[[1]][2:5], collapse = "."))
  
  # Analyze survival accross trna expression
  pdf(sprintf("plots/survival/%s_survival_RAA.pdf",type))
  for (trna in filtered_trnas){
    survInfo$trna_bin = binarize(as.numeric(tempdata[trna,]),
                                 substr(colnames(tempdata),1,15),toupper(rownames(survInfo)),0.4)
    fit <- survfit(Surv(times, vital_status) ~ trna_bin, data = survInfo)
    pval = surv_pvalue(fit,data = survInfo)[[2]]
    if ((nrow(survInfo[survInfo$trna_bin&!is.na(survInfo$trna_bin),])>1)){
      summary[type,trna] = pval
      if (!is.na(pval)&&pval<=0.05){
        print(ggsurvplot(fit, data = survInfo, risk.table = TRUE, pval=TRUE, conf.int=TRUE,
                         title=trna))
      }
    }
  }
  dev.off()
}
summary[!is.na(summary)] <- p.adjust(summary[!is.na(summary)],method="fdr")

write.csv(summary,"results/RAA_survival_analysis.csv")