#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 10:37:57 2019

@author: xhernandez
"""

# Load modules
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def transformdata(data,transf):
    aa_idx = re.findall("i?[A-Z][a-z]{2}[A-Z]{3}",str(data.index))
    data = data.loc[aa_idx,:]
    if transf=="log":
        outdata = data.apply(np.log)
        # Remove inf values
        outdata.replace([np.inf, -np.inf], np.nan,inplace=True)
    elif transf=="arcsinh":
        outdata = data.apply(np.arcsinh)
    elif transf=="rel":
        # Compute relative data
        outdata = pd.DataFrame(columns=data.columns,index=data.index)
        aa = list(set([s[0:-3] for s in outdata.index]))
        for n in aa:
            idx = [n==s for s in [l[0:-3] for l in data.index]]
            total = data.loc[idx,:].sum()
            outdata.loc[data.index[idx],:] = data.loc[data.index[idx],:]/total
            iszero = (total==0)
            if any(iszero):
                outdata.loc[data.index[idx],iszero] = 1.0/sum(idx)
        outdata.iloc[:,:] = np.float64(outdata)
    else:
        outdata=data
        
    return outdata

#%% Load CU
codons = pd.read_csv("data/codons_table.tab", sep="\t", index_col=0)
weighted_CU = pd.read_csv("results/healthy_CU.csv")
weighted_CU.iloc[:,0] = ([codons.loc[s,"AA"]+s for s in weighted_CU.iloc[:,0]])
weighted_CU.set_index("Unnamed: 0", inplace=True)
weighted_CU.columns = [s[:-12] for s in weighted_CU.columns]
# Load trnas
trna = pd.read_csv("data/TCGAall_nomod.csv",index_col=0)
# Trnas with no copies
no_copies_trna = ['AlaGGC', 'ArgGCG', 'AspATC', 'CysACA', 'GlyACC', 'HisATG', 'LeuGAG',
                  'PheAAA', 'ProGGG', 'SerACT', 'SerGGA', 'ThrGGT', 'ValGAC']

#%% Correlate trna vs codons
# TRNA data
trna_data = transformdata(trna,"rel")
# CU data
codon_data = transformdata(weighted_CU,"rel")
anticodons = [s[:3]+codons.loc[s[-3:],"ANTICODON"] for s in codon_data.index]
# Map anticodon data to trna samples
CU = pd.DataFrame(index=anticodons,columns=trna.columns)
samples = [s[-35:-19] if (s[-1] in ["A","B"]) else s[-34:-18] for s in CU.columns]
CU.loc[:,:] = np.array([np.array(codon_data.loc[:,s].values) if (s in codon_data.columns) else [np.nan]*codon_data.shape[0] for s in samples]).transpose()

corcoef = pd.DataFrame(index=trna_data.columns,columns=trna_data.columns)
cor_anticodons = np.array(anticodons); cor_anticodons=cor_anticodons[np.array([not s in no_copies_trna for s in cor_anticodons])]
for sample in trna_data.columns:
    corDf = pd.concat([trna_data.loc[cor_anticodons,sample],CU.loc[cor_anticodons,:]],axis=1)
    corcoef.loc[sample,:] = corDf.corr(method="spearman").iloc[0,1:]

corcoef.to_csv("results/corcoef_relcodons.csv")

#%% Corcoef data
corcoef = pd.read_csv("results/corcoef_relcodons.csv",index_col=0)
