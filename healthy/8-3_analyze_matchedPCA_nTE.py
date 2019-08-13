# -*- coding: utf-8 -*-

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
trna = pd.read_csv("results/matched_AAcTE.csv",index_col=0).transpose()

mrnaseq = pd.read_csv("data/healthy_mrnaseq.csv")
mrnaseq.set_index('Unnamed: 0',inplace=True)
mrnaseq.columns = [s[:-12] for s in mrnaseq.columns]

#%% Dimension reduction
data = transformdata(trna.dropna(axis = 1, how = 'any'),"") #Remove NaN

## PCA
pca = PCA(n_components=2)
x = StandardScaler().fit_transform(data.transpose())
principalComponents = pca.fit_transform(x)
expl_var = pca.explained_variance_ratio_
features = pd.DataFrame(pca.components_.transpose(), columns=["PCA1","PCA2"], index = data.index)

## Plot
principalDf = pd.DataFrame(principalComponents,columns=["PCA1","PCA2"],index=data.columns)
principalDf["source"] = [re.findall("[A-Za-z0-9_+]+(?=.)",s)[0] for s in data.columns]

#%% Calculate mitotic index
genes = ["MKI67"] # ki67
mitindex = pd.DataFrame(columns=["mitindex"]+genes,index=data.columns)
samples = [s[-35:-19] if (s[-1] in ["A","B"]) else s[-34:-18] for s in mitindex.index]
samples = [re.sub("\\.","-",s) for s in samples]
mitindex.loc[:,genes] = [list(mrnaseq.loc[genes,s].values) if (s in mrnaseq.columns) else [np.nan]*len(genes) for s in samples]
mitindex.loc[:,"mitindex"] = mitindex.loc[:,genes].mean(1)

#%% Plot cancer types ordered by mitotic index
fig = plt.figure(figsize = (12,10))
ax = plt.subplot() 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
ax.set_title('Dimension Reduction', fontsize = 20)

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
    if label in ["COAD", "READ", "GBM"]:#sum(idx)>1:
        for lab in notnaDf.index[idx]:
            ax.annotate(lab,[notnaDf.loc[lab, 'PCA1'],notnaDf.loc[lab, 'PCA2']],fontsize=8)

fig.legend(handles,labels,loc="center right")
ax.grid()

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
    if label in ["COAD", "READ", "GBM"]:#sum(idx)>1:
        for lab in principalDf.index[idx]:
            ax.annotate(lab,[principalDf.loc[lab, 'PCA1'],principalDf.loc[lab, 'PCA2']],fontsize=8)
ax.text(0.05,0.95,str(r"$R_{PC1-spearman} = %0.3f$; $R_{PC2-spearman} = %0.3f$" % (corcoef[0],corcoef[1])),fontsize=12,transform=ax.transAxes)
ax.grid()

#%% How does PCA1 correlate with other genes?
genes = list(mrnaseq.index)
genexp = pd.DataFrame(columns=genes,index=data.columns)
samples = [s[-35:-19] if (s[-1] in ["A","B"]) else s[-34:-18] for s in genexp.index]
samples = [re.sub("\\.","-",s) for s in samples]
genexp.loc[:,genes] = [list(mrnaseq.loc[genes,s].values) if (s in mrnaseq.columns) else [np.nan]*len(genes) for s in samples]

# Correlation
genexp["PCA1"] = principalDf.loc[:,"PCA1"]
corcoef = genexp.corr(method="spearman").iloc[:-1,-1]
#corcoef = pd.read_csv("results/spearman_correlation_PC1.csv",header=None,index_col=0)
# Randomization test 
pval = sum(corcoef.values>corcoef.loc["MKI67"])/(corcoef.shape[0]-1)