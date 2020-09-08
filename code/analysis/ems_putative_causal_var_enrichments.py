import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')

import time as tm
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

#read the enrichment dataframe
fn = "gs://qingbowang/tmp/ukb_100k_allscores.tsv"
with hl.hadoop_open(fn, 'w') as f:
    df1.to_csv(f, sep="\t")
fn = "gs://qingbowang/tmp/gtex_100k_allscores.tsv"
with hl.hadoop_open(fn, 'w') as f:
    df2.to_csv(f, sep="\t")
fn = "gs://qingbowang/tmp/bbj_100k_allscores.tsv"
with hl.hadoop_open(fn, 'w') as f:
    df3.to_csv(f, sep="\t")
fn = "gs://qingbowang/tmp/gev_100k_allscores.tsv"
with hl.hadoop_open(fn, 'w') as f:
    df4.to_csv(f, sep="\t")
    
df3["putative_causal"] = df3.max_pip>0.9
df4["putative_causal"] = df4.pip_dap>0.9

F1 = sum(df1.putative_causal)/df1.shape[0]
F2 = sum(df2.putative_causal)/df2.shape[0]
F3 = sum(df3.putative_causal)/df3.shape[0]
F4 = sum(df4.putative_causal)/df4.shape[0]

#plot GTEx and UKB enrichments:
colors = ["tab:green", "tab:orange", "tab:blue", "tab:pink", "tab:olive", "tab:cyan", "tab:brown"]

df1["abs_gerp"] = abs(df1.gerp)
df1["tss"] = 1/df1.min_tss

df2["abs_gerp"] = abs(df2.gerp)
df2["tss"] = 1/abs(df2.tss_distance)

plt.rcParams.update({'font.size': 18})
fig, (ax0, ax1) = plt.subplots(2,1, figsize=(12,6))
i = 0
#plot eqtls:

for s in ["ems","tss","mean_e_val","cadd","ncer","abs_gerp","Non-Coding_Score"]:
    df2["x_bin"] = abs(df2[s].astype(float)).rank(pct=True) // 0.01
    #df2["x_bin"] = pd.qcut(abs(df2[s].astype(float)), 100, labels=False)
    y_bin = df2[~np.isnan(df2.x_bin)].groupby("x_bin").putative_causal.agg(['mean','sem']) / F2
    #impute in case having NAs: just happens once in gerp
    setd = np.setdiff1d(np.arange(100), y_bin.index)
    if len(setd)>0:
        for b in setd:
            y_bin.loc[b,:] = np.nan  
    #ax0.errorbar((np.arange(100)-50)/120 + i, y_bin["mean"], yerr=y_bin["sem"], fmt="o", color=colors[i])
    ax0.errorbar((y_bin.index-50)/120 + i, y_bin["mean"], yerr=y_bin["sem"], fmt="o", color=colors[i])
    i += 1
ax0.set_xticklabels([])
    
#plot comp trait:
i = 0
for s in ["ems","tss","mean_e_val","cadd","ncer", "abs_gerp","Non-Coding_Score"]:
    df1["x_bin"] = abs(df1[s].astype(float)).rank(pct=True) // 0.01
    #df1["x_bin"] = pd.qcut(abs(df1[s].astype(float)), 100, labels=False)
    y_bin = df1[~np.isnan(df1.x_bin)].groupby("x_bin").putative_causal.agg(['mean','sem']) / F1
    #ax1.errorbar((np.arange(100)-50)/120 + i, y_bin["mean"], yerr=y_bin["sem"], fmt="o", color=colors[i])
    ax1.errorbar((y_bin.index-50)/120 + i, y_bin["mean"], yerr=y_bin["sem"], fmt="o", color=colors[i])
    i += 1
#ax0.set_ylabel("dataset: GTEx\ntrait: whole blood eQTL", size=18)    
#ax1.set_ylabel("dataset: UKBB\ntrait: Hematopoietic", size=18)
#ax0.text(2.6,15.5,"dataset: GTEx, trait: whole blood eQTL", color="tab:green")
#ax1.text(2.8,11, "dataset: UKBB, trait: hematopoietic", color="tab:red")

ax0.set_ylabel("GTEx\nenrichment", fontsize=23)
ax1.set_ylabel("UKBB\nenrichment", fontsize=23)
plt.xticks(np.arange(7), ["EMS", "TSS", "DeepSEA","CADD", "ncER", "GERP","Fathmm"], size=23, rotation=15)
plt.xlabel("Score (percentile)", size=24)
#fig.text(0, 0.55, 'putative causal variant enrichment', va='center', rotation='vertical', fontsize=16)
fig.tight_layout()
fn = "gs://qingbowang/ems_v1_test/fig2/fig2_enrichments.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()

#plot BBJ and Geuvadis
df3["abs_gerp"] = abs(df3.gerp)
df3["tss"] = 1/df3.min_abs_tss_dist
df4["abs_gerp"] = abs(df4.gerp)
df4["tss"] = 1/abs(df4.tss_distance)

plt.rcParams.update({'font.size': 18})
fig, (ax0, ax1) = plt.subplots(2,1, figsize=(12,6))
i = 0
#plot eqtls:
for s in ["ems","tss","mean_e_val","cadd","ncer", "abs_gerp","Non-Coding_Score"]:
    df4["x_bin"] = abs(df4[s]).rank(pct=True) // 0.01
    y_bin = df4[~np.isnan(df4.x_bin)].groupby("x_bin").putative_causal.agg(['mean','sem']) / F4
    #impute in case having NAs: just happens once in gerp
    setd = np.setdiff1d(np.arange(100), y_bin.index)
    if len(setd)>0:
        for b in setd:
            y_bin.loc[b,:] = np.nan  
    ax0.errorbar((y_bin.index-50)/120 + i, y_bin["mean"], yerr=y_bin["sem"], fmt="o", color=colors[i])
    i += 1
ax0.set_xticklabels([])
    
#plot comp trait:
i = 0
for s in ["max_ems","tss","mean_e_val","cadd","ncer", "abs_gerp","Non-Coding_Score"]:
    df3["x_bin"] = abs(df3[s]).rank(pct=True) // 0.01
    y_bin = df3[~np.isnan(df3.x_bin)].groupby("x_bin").putative_causal.agg(['mean','sem']) / F3
    ax1.errorbar((y_bin.index-50)/120 + i, y_bin["mean"], yerr=y_bin["sem"], fmt="o", color=colors[i])
    i += 1
#ax0.text(3.01,21,"dataset: Geuvadis, trait: LCL eQTL", color="tab:green")
#ax1.text(3.01,15, "dataset: BBJ, trait: hematopoietic", color="tab:red")
ax0.set_ylabel("Geuvadis\nenrichment", fontsize=23)
ax1.set_ylabel("BBJ\nenrichment", fontsize=23)

plt.xticks(np.arange(7), ["EMS", "TSS", "DeepSEA","CADD", "ncER", "GERP","Fathmm"], size=22, rotation=15)
plt.xlabel("Score (percentile)", size=22)
#fig.text(0, 0.55, 'putative causal variant enrichment', va='center', rotation='vertical', fontsize=16)
fig.tight_layout()
fn = "gs://qingbowang/ems_v1_test/fig2/fig2_enrichment_replication.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()
