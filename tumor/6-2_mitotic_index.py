#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 15:44:29 2019

@author: xhernandez
"""

# Load modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
#from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def transformdata(data,transf):
    data = data.iloc[0:66,:]
    if transf=="log":
        outdata = data.apply(np.log)
        # Remove inf values
        outdata.replace([np.inf, -np.inf], np.nan,inplace=True)
    elif transf=="arcsinh":
        outdata = data.apply(np.arcsinh)
    elif transf=="sqrt":
        outdata = data.apply(np.sqrt)
    elif transf=="rel":
        data = data.iloc[0:61,:]
        # Compute relative data
        outdata = pd.DataFrame(columns=data.columns,index=data.index[0:61])
        aa = list(set([s[0:3] for s in outdata.index]))
        for n in aa:
            idx = [n in s for s in [l[0:3] for l in data.index]]
            total = data.loc[idx,:].sum()
            outdata.loc[data.index[idx],:] = data.loc[data.index[idx],:]/total
            iszero = (total==0)
            if any(iszero):
                outdata.loc[data.index[idx],iszero] = 1.0/sum(idx)
        outdata.iloc[:,:] = np.float64(outdata)
    else:
        outdata=data
        
    return outdata

def get_expression(abbr):
    if ";" in abbr:
        names = abbr.split(";")
        output = pd.DataFrame()
        for n in names:
            new = pd.read_csv(str("/home/xhernandez/Downloads/TCGA-mRNAseq/20160128-%s-RNAseq2GeneNorm.txt" % n), 
                             index_col=0, sep= "\t",na_values="normalized_count")
            rownames= np.array([s.split("|")[0] for s in new.index])
            idx = np.array([s not in ["gene_id", "?"] for s in rownames])
            new = new.loc[idx,:].set_index(rownames[idx])
            new.rename(lambda x: x[:16],axis=1,inplace=True)
            output = pd.concat([output,new], axis=1)
    else:
        output = pd.read_csv(str("/home/xhernandez/Downloads/TCGA-mRNAseq/20160128-%s-RNAseq2GeneNorm.txt" % abbr), 
                             index_col=0, sep= "\t",na_values="normalized_count")
        rownames= np.array([s.split("|")[0] for s in output.index])
        idx = np.array([s not in ["gene_id", "?"] for s in rownames])
        output = output.loc[idx,:].set_index(rownames[idx])
        output.rename(lambda x: x[:16],axis=1,inplace=True)
        
    return output

#%% Load data
cancer_types = {"BRCA":"BRCA","PRAD":"PRAD","kidney":"KICH;KIRP;KIRC","lung":"LUAD;LUSC","HNSC":"HNSC","uterus":"UCEC;CESC",
"liver":"LIHC;CHOL","THCA":"THCA","colorectal":"COAD;READ","ESCA":"ESCA","STAD":"STAD","BLCA":"BLCA","PAAD":"PAAD","THYM":"THYM",
"SKCM":"SKCM","PCPG":"PCPG"}

trna = pd.read_csv("results/AAcTE_CUgenomic_sqrt.csv",index_col=0)

#%% Load expression data
mrnaseq = None
for catype in cancer_types.keys():
    abbr = cancer_types[catype].split(";")
    if mrnaseq is None:
        mrnaseq = get_expression(cancer_types[catype])
    else:        
        add = get_expression(cancer_types[catype])
        mrnaseq = pd.concat([mrnaseq,add], axis=1)

#%% Dimension reduction
data = transformdata(trna,"")
data = data.dropna(axis = 1, how = 'any') #Remove NaN

## PCA
pca = PCA(n_components=2)
x = StandardScaler().fit_transform(data.transpose())
principalComponents = pca.fit_transform(x)
expl_var = pca.explained_variance_ratio_
features = pd.DataFrame(pca.components_.transpose(), columns=["PCA1","PCA2"], index = data.index)


## Plot
principalDf = pd.DataFrame(principalComponents,columns=["PCA1","PCA2"],index=data.columns)
principalDf["source"] = [re.findall("[A-Za-z0-9_+]+(?=-)",s)[0] for s in data.columns]

#%% Calculate mitotic index
genes = ["MKI67"] # ki67
mitindex = pd.DataFrame(columns=["mitindex"]+genes,index=data.columns)
samples = [s[-35:-19] if (s[-1] in ["A","B"]) else s[-34:-18] for s in mitindex.index]
mitindex.loc[:,genes] = [list(mrnaseq.loc[genes,s].values) if (s in mrnaseq.columns) else [np.nan]*len(genes) for s in samples]
mitindex.loc[:,"mitindex"] = mitindex.loc[:,genes].mean(1)

#%% Plot mitotic index
fig = plt.figure(figsize = (15,10))
ax = plt.subplot() 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
ax.set_title('Dimension Reduction', fontsize = 20)

plt.scatter(principalDf.loc[:, 'PCA1'], principalDf.loc[:, 'PCA2'], c = mitindex.loc[:,"mitindex"].apply(np.log), s = 50, alpha=0.5, cmap="jet")
cb = plt.colorbar()
#cb.set_label('Mitotic Index')
corDf1 = pd.DataFrame([principalDf.loc[:,"PCA1"],mitindex.loc[:,"mitindex"]]).transpose()
corDf2 = pd.DataFrame([principalDf.loc[:,"PCA2"],mitindex.loc[:,"mitindex"]]).transpose()
corcoef = [corDf1.corr(method="spearman").iloc[0,1],corDf2.corr(method="spearman").iloc[0,1]]

labels = list(set(principalDf["source"]))
for label in labels:
    idx = principalDf["source"]==label

ax.text(0.05,0.95,str(r"$R_{PC1-spearman} = %0.3f$; $R_{PC2-spearman} = %0.3f$" % (corcoef[0],corcoef[1])),fontsize=12,transform=ax.transAxes)
ax.grid()

#%% Calculate mitotic index from proteomics
cancer_types = {"colorectal":"COAD;READ"}#{"BRCA":"BRCA"}
genes = ["MKI67"] # ki67
protdata = pd.DataFrame(index=genes)
for catype in cancer_types.keys():
    add = pd.read_csv(str("/home/xhernandez/Downloads/TCGA-proteomics/TCGA_%s.csv" % catype),index_col = 0)
    protdata = pd.concat([protdata,add.loc[genes,:]],1)    
        
mitindex = pd.DataFrame(columns=["mitindex"]+genes,index=data.columns)
samples = [s[-35:-19] if (s[-1] in ["A","B"]) else s[-34:-18] for s in mitindex.index]
samples = [".".join(s.split("-")) for s in samples]
mitindex.loc[:,genes] = [list(protdata.loc[genes,s].values) if (s in protdata.columns) else [np.nan]*len(genes) for s in samples]
mitindex.loc[:,"mitindex"] = mitindex.loc[:,genes].mean(1)

#%% Plot proteomics mitotic index
fig = plt.figure(figsize = (15,10))
ax = plt.subplot() 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
ax.set_title('Dimension Reduction', fontsize = 20)

plt.scatter(principalDf.loc[:, 'PCA1'], principalDf.loc[:, 'PCA2'], c = mitindex.loc[:,"mitindex"].apply(np.log), s = 50, alpha=0.5, cmap="jet")
cb = plt.colorbar()
#cb.set_label('Mitotic Index')
corDf1 = pd.DataFrame([principalDf.loc[:,"PCA1"],mitindex.loc[:,"mitindex"]]).transpose()
corDf2 = pd.DataFrame([principalDf.loc[:,"PCA2"],mitindex.loc[:,"mitindex"]]).transpose()
corcoef = [corDf1.corr(method="spearman").iloc[0,1],corDf2.corr(method="spearman").iloc[0,1]]

labels = list(set(principalDf["source"]))
for label in labels:
    idx = principalDf["source"]==label

ax.text(0.05,0.95,str(r"$R_{PC1-spearman} = %0.3f$; $R_{PC2-spearman} = %0.3f$" % (corcoef[0],corcoef[1])),fontsize=12,transform=ax.transAxes)
ax.grid()
