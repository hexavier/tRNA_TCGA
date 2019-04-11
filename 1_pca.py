# -*- coding: utf-8 -*-

#### PCA and correlation analysis for tRNA data ####

# Load modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import re
from sklearn.decomposition import PCA
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
trna = pd.read_csv("data/HEK_nomod_norm.csv",index_col=0)
#trna = pd.read_csv("data/HEK_nomod_norm.csv",index_col=0)
no_copies_trna = np.array([False, False, True, False, True, False, False, False, False, False, False, False, True, False, 
                           True, False, False, False, False, False, False, False, False, True, False, True, False, False, 
                           False, False, False, False, False, False, True, False, False, False, False, True, False, False, 
                           True, False, True, False, False, False, True, False, False, False, False, True, False, False, 
                           False, False, False, False, True])

#%% Calculate correlation
pdf = matplotlib.backends.backend_pdf.PdfPages("plots/correlations.pdf")
cor = pd.DataFrame(columns = ["method","transf","data","correl"])
cor.loc[:,"transf"] = ["","","log","log","sqrt","sqrt","rel","rel"]
cor.loc[:,"method"] = ["pearson","spearman","pearson","spearman","pearson","spearman","pearson","spearman"]
cmap = plt.cm.get_cmap('bwr') # set color map
normcol = plt.Normalize(-1,1)
for n in range(cor.shape[0]):
    cor.iloc[n].data = transformdata(trna,cor.iloc[n,1])
    if cor.iloc[n,1] == "rel":
        # Remove 0 and 1
        unique_aa = [(sum([s[0:3] in n for n in cor.iloc[n,2].index[~no_copies_trna]])>1) for s in cor.iloc[n,2].index]
        cor.iloc[n].data = cor.iloc[n,2].loc[unique_aa,:]
    cor.iloc[n].correl = cor.iloc[n,2].corr(method=cor.iloc[n,0])
    
    axes = pd.plotting.scatter_matrix(cor.iloc[n,2], figsize=(10, 10))
    for i, j in zip(*plt.np.triu_indices_from(axes, k=1)):
        R = cor.iloc[n,3].iloc[i,j]
        axes[i, j].annotate("%.3f" % R, (0.8, 0.8), xycoords='axes fraction', ha='center', va='center')
        axes[i, j].set_facecolor(cmap(normcol(R)))
    for ax in axes.ravel():
        ax.set_xlabel(ax.get_xlabel(), rotation = 45, horizontalalignment="right")
        ax.set_ylabel(ax.get_ylabel(), rotation = 45, horizontalalignment="right")
    fig = plt.gcf()
    fig.suptitle(str("%s transormation - %s" % (cor.iloc[n,1],cor.iloc[n,0])))
    pdf.savefig(fig, bbox_inches="tight" , pad_inches=0.2)
pdf.close()

#%% Dimension reduction
data = transformdata(trna,"rel")
data = data.dropna(axis = 0, how = 'any') #Remove NaN

## PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(data.transpose())
expl_var = pca.explained_variance_ratio_
features = pd.DataFrame(pca.components_.transpose(), columns=["PCA1","PCA2"], index = data.index)

## tSNE
#principalComponents = TSNE(n_components=2,n_iter=5000,perplexity=30).fit_transform(data.transpose())
#expl_var = [0,0]

## Plot
principalDf = pd.DataFrame(principalComponents,columns=["PCA1","PCA2"],index=data.columns)
principalDf["source"] = [re.findall("[A-Za-z0-9_+]+(?=-)",s)[0] for s in data.columns]
#principalDf["source"] = ["small-RNAseq","small-RNAseq","small-RNAseq","small-RNAseq","tRNA-Seq",
#           "tRNA-Seq","tRNA-Seq","copy number","other","other","other","other","other","other","other","other"]
#principalDf["source"] = ["small-RNAseq","small-RNAseq","small-RNAseq","small-RNAseq","precursors","precursors","precursors","precursors"]


# Plot pca
fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
ax.set_title('Dimension Reduction', fontsize = 20)

labels = list(set(principalDf["source"]))
cmap = plt.cm.get_cmap('jet') # set color map
colors = cmap(np.arange(len(labels))/(len(labels)-1))
handles = [mlines.Line2D([], [], marker='o', color=c) for c in colors]

for label, color in zip(labels,colors):
    idx = principalDf["source"]==label
    ax.scatter(principalDf.loc[idx, 'PCA1'], principalDf.loc[idx, 'PCA2'], c = color, s = 50, alpha=0.5)
    if sum(idx)>1:#label in ["COAD", "READ", "GBM"]:
        for lab in principalDf.index[idx]:
            ax.annotate(lab,[principalDf.loc[lab, 'PCA1'],principalDf.loc[lab, 'PCA2']],fontsize=8)
fig.legend(handles,labels,loc="center right")
ax.grid()
