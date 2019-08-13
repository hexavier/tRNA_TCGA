#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 10:35:09 2019

@author: xhernandez
"""

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

#%% Load data
trna = pd.read_csv("data/TCGAall_nomod.csv",index_col=0)
#trna = pd.read_csv("data/HEK_nomod_norm.csv",index_col=0)
mrnaseq = pd.read_csv("data/healthy_mrnaseq.csv")
mrnaseq.set_index('Unnamed: 0',inplace=True)
mrnaseq.columns = [s[:-12] for s in mrnaseq.columns]

#%% Calculate mitotic index
prolif_genes = ["MKI67"]
mitindex = pd.DataFrame(columns=["mitindex"]+prolif_genes,index=trna.columns)
samples = [s[-35:-19] if (s[-1] in ["A","B"]) else s[-34:-18] for s in mitindex.index]
mitindex.loc[:,prolif_genes] = [list(mrnaseq.loc[prolif_genes,s].values) if (s in mrnaseq.columns) else [np.nan]*len(prolif_genes) for s in samples]
mitindex.loc[:,"mitindex"] = mitindex.loc[:,prolif_genes].mean(1)

#%% Dimension reduction
data = transformdata(trna,"rel")
data = data.dropna(axis = 0, how = 'any') #Remove NaN
targets=[re.findall("[A-Za-z0-9_+]+(?=-)",s)[0] for s in data.columns]

## Discriminant analysis
disc = LinearDiscriminantAnalysis(n_components=2)
principalComponents = disc.fit_transform(data.transpose(),targets)
expl_var = disc.explained_variance_ratio_
coef = pd.DataFrame(disc.coef_.transpose(), columns=disc.classes_, index = data.index)

## Plot
principalDf = pd.DataFrame(principalComponents,columns=["PCA1","PCA2"],index=data.columns)
principalDf["source"] = targets

#%% Plot cancer types ordered by mitotic index
fig = plt.figure(figsize = (12,10))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
ax.set_title('Discriminant Analysis', fontsize = 20)

notnaDf = principalDf.loc[mitindex.loc[:,"mitindex"].notna(),:]
mitnotna = mitindex.dropna(0)
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
    if label in ["COAD", "READ", "GBM"]:#sum(idx)>1:
        for lab in notnaDf.index[idx]:
            ax.annotate(lab,[notnaDf.loc[lab, 'PCA1'],notnaDf.loc[lab, 'PCA2']],fontsize=8)

fig.legend(handles,labels,loc="center right")
ax.grid()

#%% Plot cancer types
fig = plt.figure(figsize = (12,10))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
ax.set_title('Discriminant Analysis', fontsize = 20)

#labels = list(set(principalDf["source"]))
cmap = plt.cm.get_cmap('jet') # set color map
colors = cmap(np.arange(len(labels))/(len(labels)-1))
handles = [mlines.Line2D([], [], marker='o', color=c, alpha=0.5) for c in colors]

for label, color in zip(labels,colors):
    idx = principalDf["source"]==label
    ax.scatter(principalDf.loc[idx, 'PCA1'], principalDf.loc[idx, 'PCA2'], c = color, s = 50, alpha=0.5)
    if label in ["COAD", "READ", "GBM"]:#sum(idx)>1:
        for lab in principalDf.index[idx]:
            ax.annotate(lab,[principalDf.loc[lab, 'PCA1'],principalDf.loc[lab, 'PCA2']],fontsize=8)
fig.legend(handles,labels,loc="center right")
ax.grid()

#%% Plot mitotic index
fig = plt.figure(figsize = (15,10))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
ax.set_title('Discriminant Analysis', fontsize = 20)

plt.scatter(principalDf.loc[:, 'PCA1'], principalDf.loc[:, 'PCA2'], c = mitindex.loc[:,"mitindex"].apply(np.log), s = 50, alpha=0.5, cmap="jet")
cb = plt.colorbar()
cb.set_label('Mitotic Index')

labels = list(set(principalDf["source"]))
for label in labels:
    idx = principalDf["source"]==label
    if label in ["COAD", "READ", "GBM"]:#sum(idx)>1:
        for lab in principalDf.index[idx]:
            ax.annotate(lab,[principalDf.loc[lab, 'PCA1'],principalDf.loc[lab, 'PCA2']],fontsize=8)
ax.grid()
