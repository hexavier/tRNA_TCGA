#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue June 04 15:44:29 2019

@author: xhernandez
"""

# Load modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import re
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import matplotlib.lines as mlines

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

trnaH = pd.read_csv("data/TCGAhealthy_nomod.csv",index_col=0)
trnaT = pd.read_csv("data/TCGAtumor_nomod.csv",index_col=0)
trna = pd.concat([trnaH,trnaT],axis=1)

#%%
pdf = matplotlib.backends.backend_pdf.PdfPages("plots/discriminant_analyses_rel.pdf")
weights = pd.DataFrame(columns=cancer_types.keys(), index=trna.index)

for catype in cancer_types.keys():
    abbr = cancer_types[catype].split(";")
    # Arrange tumor type data
    caidx = np.array([s in abbr for s in [re.findall("[A-Za-z0-9_+]+(?=-)",s)[0] for s in trna.columns]])
    data = transformdata(trna.loc[:,caidx],"rel")
    data = data.dropna(axis = 0, how = 'any') #Remove NaN
    
    mrnaseq = get_expression(cancer_types[catype])

    #%% Calculate mitotic index
    genes = ["MKI67"] # ki67
    mitindex = pd.DataFrame(columns=["mitindex"]+genes,index=data.columns)
    samples = [s[-35:-19] if (s[-1] in ["A","B"]) else s[-34:-18] for s in mitindex.index]
    mitindex.loc[:,genes] = [list(mrnaseq.loc[genes,s].values) if (s in mrnaseq.columns) else [np.nan]*len(genes) for s in samples]
    mitindex.loc[:,"mitindex"] = mitindex.loc[:,genes].mean(1)
    
    #%% Dimension reduction
    ## Discriminant analysis
    tumcode = [re.findall("(?<=TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-)[0-9]+",s)[0] for s in data.columns]
    targets = ["CA" if s=="01" else "HE" for s in tumcode]
        
    ## Discriminant analysis
    disc = LinearDiscriminantAnalysis(n_components=2, store_covariance=True)
    principalComponents = disc.fit_transform(data.transpose(),targets)
    expl_var = disc.explained_variance_ratio_
    normcoef = disc.coef_*disc.covariance_.diagonal()
    coef = pd.DataFrame(normcoef.transpose(), index = data.index)
    # Record normalized coefficients
    weights.loc[coef.index,catype] = coef.values.reshape((coef.shape[0],))
    
    # Add an uninformative second component if not generated
    if principalComponents.shape[1]==1:
        principalComponents = np.append(principalComponents,np.random.uniform(size=(principalComponents.shape[0],1)),axis=1)
        expl_var = np.append(expl_var, 0)

    ## Plot
    principalDf = pd.DataFrame(principalComponents,columns=["PCA1","PCA2"],index=data.columns)
    principalDf["source"] = targets

    #%% Plot cancer types ordered by mitotic index
    fig = plt.figure(figsize = (12,10))
    ax = plt.subplot() 
    ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
    ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
    ax.set_title(str('Dimension Reduction %s' % catype), fontsize = 20)
    
    notnaDf = principalDf.loc[mitindex.loc[:,"mitindex"].notna(),:]
    mitnotna = mitindex.dropna(0,how="all")
    labels = np.array(list(set(notnaDf["source"])))
    
    means=[]
    for label in labels:
        idx = notnaDf["source"]==label
        means.append(mitnotna.loc[idx,"mitindex"].apply(np.log).mean(0))
    meanssort = np.array(means); meanssort.sort()
    order = np.array([[s==n for n in means].index(True) for s in meanssort])
    labels = labels[order]
    
    cmap = plt.cm.get_cmap('jet') # set color map
    normcol = plt.Normalize(min(meanssort),max(meanssort))
    colors = cmap(normcol(meanssort))#cmap(np.arange(len(labels))/(len(labels)-1))
    handles = [mlines.Line2D([], [], marker='o', color=c, alpha=0.5) for c in colors]
    
    for label, color in zip(labels,colors):
        idx = notnaDf["source"]==label
        ax.scatter(notnaDf.loc[idx, 'PCA1'], notnaDf.loc[idx, 'PCA2'], c = color, s = 50, alpha=0.5)
    
    fig.legend(handles,labels,loc="center right")
    ax.grid()
    
    pdf.savefig(fig, bbox_inches="tight" , pad_inches=0.2)
    
pdf.close()
weights.to_csv("results/discriminant_analyses_rel.csv")