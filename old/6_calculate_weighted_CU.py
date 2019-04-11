#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 17:04:05 2019

@author: xhernandez
"""

# Load modules
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def transformdata(data,transf):
    aa_idx = re.findall("[A-Z][a-z]{2}[A-Z]{3}",str(data.index))
    data = data.loc[aa_idx,:]
    if transf=="log":
        outdata = data.apply(np.log)
        # Remove inf values
        outdata.replace([np.inf, -np.inf], np.nan,inplace=True)
    elif transf=="arcsinh":
        outdata = data.apply(np.arcsinh)
    elif transf=="rel":
        data = data.iloc[0:64,:]
        # Compute relative data
        outdata = pd.DataFrame(columns=data.columns,index=data.index[0:64])
        aa = list(set([s[0:-3] for s in outdata.index]))
        for n in aa:
            idx = [n in s for s in [l[0:-3] for l in data.index]]
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
# Genomic codon usage https://hive.biochemistry.gwu.edu/dna.cgi?cmd=refseq_processor&id=569942
codus = pd.read_csv("data/human_CU_refseq.tsv",sep="\t",index_col="Protein ID")

# Isoform gene expression
isofseq = pd.read_csv("data/healthy_isoform_seq.csv")
isofseq.loc[:,"Unnamed: 0"] = [s[:-2] for s in isofseq.loc[:,"Unnamed: 0"]]
isofseq.set_index("Unnamed: 0",inplace=True)

# Mapping
tcga2refseq = pd.read_csv("data/TCGA2RefSeq.txt",sep="\t", header=None)
# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownIsoforms.txt.gz
tcga2refseq.loc[:,0] = [s[:-2] for s in tcga2refseq.loc[:,0]]
tcga2refseq.set_index(0,inplace=True)

np2nm = pd.read_csv("data/NP2NM.txt",sep="\t")
# http://genome-euro.ucsc.edu/cgi-bin/hgTables ncbiRefLink
np2nm.loc[:,"protAcc"] = [s[:-2] if isinstance(s,str) else s for s in np2nm.loc[:,"protAcc"]]
np2nm.set_index("protAcc",inplace=True)

#%% Prepare CU data
# Map ids
tcga_ids = np.array([s for s in isofseq.index if s in tcga2refseq.index])
refseq_ids = np.array([tcga2refseq.loc[s,1] for s in isofseq.index if s in tcga2refseq.index])
# Detect whether refseq ids are in codon usage data, and remove if not
codusindex_main = np.array([np2nm.loc[n[:-2],"#id"][:-2] if n[:-2] in np2nm.index else n[:-2] for n in codus.index])

# Keep only events in common between tcga and codon usage
keep = np.array([s in codusindex_main for s in refseq_ids])
tcga_ids = tcga_ids[keep]
refseq_ids = refseq_ids[keep]
# Convert booleans in codus_ids to string index. If more than 1, take mean
codus_clean = np.array([codus.loc[codusindex_main==s,:].mean(0) for s in refseq_ids])

# Clean isofseq dataset
isofseq = isofseq.loc[tcga_ids,:]
#isofseq = isofseq.apply(np.arcsinh) # transform data to avoid one gene accounting for most CU
#%% Calculate weighted CU

bins = 20
items_per_bin = int(isofseq.shape[0]/bins)
start = 0
final = items_per_bin
bybins = []
for n in range(bins):
    weighted_CU = pd.DataFrame(index=codus.columns[11:],columns=isofseq.columns)
    for i in weighted_CU.columns:
        # Sort by expression
        exprsorted = isofseq.loc[:,i].sort_values(ascending=False)
        # Calculate CU
        sample_cu = np.transpose(codus_clean[start:final,6:])*exprsorted.iloc[start:final].values
        weighted_CU.loc[:,i] = sample_cu.mean(1)
    weighted_CU.to_csv(str("results/CUbybins/healthy_CU_bin%s.csv" % n))
    start += items_per_bin
    final += items_per_bin
    if n==(bins-2):
        final = isofseq.shape[0]

#%% Load CU
codons = pd.read_csv("data/codons_table.tab", sep="\t", index_col=0)
weighted_CU = pd.read_csv("results/healthy_weightedCU.csv")
weighted_CU.iloc[:,0] = ([codons.loc[s,"AA"]+s for s in weighted_CU.iloc[:,0]])
weighted_CU.set_index("Unnamed: 0", inplace=True)

#%% Calculate mitotic index
genes = ["uc001lke","uc009yaw","uc009yav","uc001lkf"] # ki67
mitindex = pd.DataFrame(columns=["mitindex"]+genes,index=weighted_CU.columns)
samples = mitindex.index
mitindex.loc[:,genes] = [list(isofseq.loc[genes,s].values) if (s in isofseq.columns) else [np.nan]*len(genes) for s in samples]
mitindex.loc[:,"mitindex"] = mitindex.loc[:,genes].mean(1)

#%% Dimension reduction
data = transformdata(weighted_CU,"rel")
data = data.dropna(axis = 0, how = 'any') #Remove NaN

## PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(data.transpose())
expl_var = pca.explained_variance_ratio_
features = pd.DataFrame(pca.components_.transpose(), columns=["PCA1","PCA2"], index = data.index)

## Plot
principalDf = pd.DataFrame(principalComponents,columns=["PCA1","PCA2"],index=data.columns)
principalDf["source"] = [re.findall("[A-Za-z0-9_+]+(?=-)",s)[0] for s in data.columns]

#%% Plot mitotic index
fig = plt.figure(figsize = (15,10))
ax = plt.subplot() 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
ax.set_title('Dimension Reduction', fontsize = 20)

plt.scatter(principalDf.loc[:, 'PCA1'], principalDf.loc[:, 'PCA2'], c = mitindex.loc[:,"mitindex"].apply(np.log), s = 50, alpha=0.5, cmap="jet")
cb = plt.colorbar()
#cb.set_label('Mitotic Index')
corDf = pd.DataFrame([principalDf.loc[:,"PCA1"],mitindex.loc[:,"mitindex"]]).transpose()
corcoef = corDf.corr(method="spearman").iloc[0,1]

labels = list(set(principalDf["source"]))
for label in labels:
    idx = principalDf["source"]==label
    if label in ["COAD", "READ", "GBM"]:#sum(idx)>1:
        for lab in principalDf.index[idx]:
            ax.annotate(lab,[principalDf.loc[lab, 'PCA1'],principalDf.loc[lab, 'PCA2']],fontsize=8)
ax.text(0.05,0.95,str(r"$R_{spearman} = %0.3f$" % corcoef),fontsize=12,transform=ax.transAxes)
ax.grid()
