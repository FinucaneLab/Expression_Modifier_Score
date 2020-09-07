import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"


#plot the baseline features enrichments
tissue_name = "Whole_Blood"
with hl.hadoop_open("gs://qingbowang/ems_v1_test/agg_for_fig1/{0}_binfeatures_enrichments.tsv".format(tissue_name), 'r') as f:
    df = pd.read_csv(f, sep="\t")
df.index = df.pip_bin
df.columns = df.columns.str.replace("enr_", "")
dfuse = df.iloc[:,3:].T.sort_values(by=9).T
with hl.hadoop_open("gs://qingbowang/ems_v1_test/agg_for_fig1/{0}_binfeatures_errorbars.tsv".format(tissue_name), 'r') as f:
    dferr = pd.read_csv(f, sep="\t")
dferr.index = dferr.pip_bin
dferr.columns = dferr.columns.str.replace("err_", "")
dfuse = df.iloc[:,3:].T.sort_values(by=9).T
dferr = dferr.loc[dfuse.index, dfuse.columns]

with hl.hadoop_open("gs://qingbowang/ems_v1_test/agg_for_fig1/{0}_binfeatures_ns.tsv".format(tissue_name), 'r') as f:
    dfn = pd.read_csv(f, sep="\t")
ntot = dfn.n_total

x = np.arange(dfuse.shape[1])
cols = ["grey", "khaki","gold","orange", "red"]
legends = ["pip<0.01", "0.01<pip<0.1", "0.1<pip<0.5", "0.5<pip<0.9", "0.9<pip"]

xtk = pd.Series(dfuse.columns.str.split("_").str[:-1]).apply(lambda x: "_".join(x))
xtk = xtk.str.replace("UTR_3", "3'UTR").str.replace("UTR_5", "5'UTR")

plt.figure(figsize=(6,4.5))
for i in range(len(dfuse.index)):
    y = dfuse.iloc[i,:]
    yerr = dferr.iloc[i,:]
    plt.errorbar(x+i*0.075, y, yerr, color=cols[i], label=legends[i]+", n={0}".format(ntot.iloc[i]), fmt="o")
plt.xticks(x, xtk, rotation=60)
plt.xlabel("functional category")
plt.ylabel("enrichment")
plt.legend(title="Whole Blood pip bin", loc="upper left")
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/fig1_baseline.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show() 

#Basenji feature:
bas_name = "CAGE:Neutrophils,"
tissue_name = "Whole_Blood"
bas_id = "CNhs10862" #corresponding to cage neutrophils
with hl.hadoop_open("gs://qingbowang/ems_v1_test/agg_for_fig1/{0}_{1}_basenji_agg_forplot.tsv".format(tissue_name, bas_id), 'r') as f:
    df = pd.read_csv(f, sep="\t", index_col=0)
nsum = pd.DataFrame(df.groupby("pip_bin").n.sum())
nsum["frac"] = nsum.n / nsum.n.sum()
df = df[~df.index.isna()]
dfall = []
for i in range(9+1):
                dfsub = df[df.index==i]
                nall = dfsub.n.sum()
                dfsub["frac"] = dfsub.n / nall
                dfsub["frac_err"] = np.sqrt(dfsub.frac*(1-dfsub.frac) / dfsub.n.sum())
                if dfsub.shape[0]<5: #if some are missing, impute as 0
                    ind = dfsub.index[0]
                    missing_pip_bins = np.setdiff1d([-1,0,1,5,9], dfsub.pip_bin)
                    for missing_pip_bin in missing_pip_bins:
                        dfsub.loc[missing_pip_bin+100] = 0 #dummy for loc 
                        dfsub.loc[missing_pip_bin+100,"pip_bin"] = missing_pip_bin
                    dfsub.sort_values(by="pip_bin", inplace=True)
                    dfsub.index = [ind]*dfsub.shape[0]
                    print ("imputing {0}, {1}, {2}".format(tissue_name, bas_id, i))
                else:
                    dfsub["frac_enr"]  = np.array(dfsub.frac) / np.array(nsum.frac)
                    dfsub["frac_enr_err"]  = np.array(dfsub.frac_err) / np.array(nsum.frac)
                dfall.append(dfsub)
dfall = pd.concat(dfall)
plt.figure(figsize=(6, 4.5))
df0 = dfall[dfall.pip_bin==-1]
plt.errorbar(df0.index, df0.frac_enr, yerr = df0.frac_enr_err, fmt="o", color="grey", label="PIP<0.01")
df001 = dfall[dfall.pip_bin==0]
plt.errorbar(df001.index+0.1, df001.frac_enr, yerr = df001.frac_enr_err, fmt="o", color="khaki", label="0.01<PIP<0.1")
df01 = dfall[dfall.pip_bin==1]
plt.errorbar(df01.index+0.2, df01.frac_enr, yerr = df01.frac_enr_err, fmt="o", color="gold", label="0.1<PIP<0.5")
df05 = dfall[dfall.pip_bin==5]
plt.errorbar(df05.index+0.3, df05.frac_enr, yerr = df05.frac_enr_err, fmt="o", color="orange", label="0.5<PIP<0.9")
df09 = dfall[dfall.pip_bin==9]
plt.errorbar(df09.index+0.4, df09.frac_enr, yerr = df09.frac_enr_err, fmt="o", color="red", label="0.9<PIP")
plt.legend(title = "Whole Blood PIP bin", loc="upper left")#, fontsize=12)
plt.xlabel("{0} Basenji score bin".format(bas_name))
plt.ylabel("Enrichment")
plt.xticks([0,9], ["<1", "100<"], rotation=30)
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/fig1_basenji_wb.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()
        
#TSS:
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


with hl.hadoop_open("gs://qingbowang/ems_v1_test/agg_for_fig1/{0}_tss_dist_agg.tsv".format(tissue_name), 'r') as f:
    df = pd.read_csv(f, sep="\t", index_col=0)
nsum = pd.DataFrame(df.groupby("pip_bin").n.sum())
nsum["frac"] = nsum.n / nsum.n.sum()
dfall = []
for i in df.index.unique():
                dfsub = df[df.index==i]
                nall = dfsub.n.sum()
                dfsub["frac"] = dfsub.n / nall
                dfsub["frac_err"] = np.sqrt(dfsub.frac*(1-dfsub.frac) / dfsub.n.sum())
                dfsub["frac_enr"]  = np.array(dfsub.frac) / np.array(nsum.frac)
                dfsub["frac_enr_err"]  = np.array(dfsub.frac_err) / np.array(nsum.frac)
                dfall.append(dfsub)
dfall = pd.concat(dfall)

xtk = ["<100", "$10^3$", "$10^4$","$10^5$","$10^6$<"]
plt.figure(figsize=(6, 4.5))
df0 = dfall[dfall.pip_bin==-1]
plt.errorbar(df0.index, df0.frac_enr, yerr = df0.frac_enr_err, fmt="o", color="grey", label="PIP<0.01")
df001 = dfall[dfall.pip_bin==0]
plt.errorbar(df001.index+0.05, df001.frac_enr, yerr = df001.frac_enr_err, fmt="o", color="khaki", label="0.01<PIP<0.1")
df01 = dfall[dfall.pip_bin==1]
plt.errorbar(df01.index+0.1, df01.frac_enr, yerr = df01.frac_enr_err, fmt="o", color="gold", label="0.1<PIP<0.5")
df05 = dfall[dfall.pip_bin==5]
plt.errorbar(df05.index+0.15, df05.frac_enr, yerr = df05.frac_enr_err, fmt="o", color="orange", label="0.5<PIP<0.9")
df09 = dfall[dfall.pip_bin==9]
plt.errorbar(df09.index+0.2, df09.frac_enr, yerr = df09.frac_enr_err, fmt="o", color="red", label="0.9<PIP")
plt.legend(title = "Whole Blood PIP bin", loc="upper right")#, fontsize=12)
plt.xlabel("distance to TSS bin")
plt.ylabel("Enrichment")
plt.xticks([2,3,4,5,6],xtk, rotation=30)
plt.tight_layout()

fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/fig1_tssdist_wb.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)

#fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/tssdist_vs_wb_pip.pdf"
#with hl.hadoop_open(fn, 'wb') as f:
#    plt.savefig(f, dpi=400)

plt.show()


#EMS itself:
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/{0}_ems_agg.tsv".format(tissue_name)
with hl.hadoop_open(fn, 'r') as f:
    df = pd.read_csv(f, sep="\t")
df = df[~df.pip_bin.isna()]    
df = df[~df.ems_bin.isna()]    
nsum = pd.DataFrame(df.groupby("pip_bin").n.sum())
nsum["frac"] = nsum.n / nsum.n.sum()
dfall = []
for i in df.ems_bin.unique():
                dfsub = df[df.ems_bin==i]
                nall = dfsub.n.sum()
                dfsub["frac"] = dfsub.n / nall
                dfsub["frac_err"] = np.sqrt(dfsub.frac*(1-dfsub.frac) / dfsub.n.sum())
                dfsub["frac_enr"]  = np.array(dfsub.frac) / np.array(nsum.frac)
                dfsub["frac_enr_err"]  = np.array(dfsub.frac_err) / np.array(nsum.frac)
                dfall.append(dfsub)
dfall = pd.concat(dfall)
nsum2 = df.groupby("ems_bin").n.sum()
xtk = ["<0.1\nn={0}".format(nsum2[-1]),"(0.1,1]\nn={0}".format(nsum2[0]),
       "[1,10)\nn={0}".format(nsum2[1]),"[10,100)\nn={0}".format(nsum2[2]),
       "[100,1000)\nn={0}".format(nsum2[3]),"1000<\nn={0}".format(nsum2[4])]
plt.figure(figsize=(6, 4.5))
df0 = dfall[dfall.pip_bin==-1]
plt.errorbar(np.arange(6), df0.frac_enr, yerr = df0.frac_enr_err, fmt="o", color="grey", label="PIP<0.01")
df001 = dfall[dfall.pip_bin==0]
plt.errorbar(np.arange(6)+0.1, df001.frac_enr, yerr = df001.frac_enr_err, fmt="o", color="khaki", label="0.01<PIP<0.1")
df01 = dfall[dfall.pip_bin==1]
plt.errorbar(np.arange(6)+0.2, df01.frac_enr, yerr = df01.frac_enr_err, fmt="o", color="gold", label="0.1<PIP<0.5")
df05 = dfall[dfall.pip_bin==5]
plt.errorbar(np.arange(6)+0.3, df05.frac_enr, yerr = df05.frac_enr_err, fmt="o", color="orange", label="0.5<PIP<0.9")
df09 = dfall[dfall.pip_bin==9]
plt.errorbar(np.arange(6)+0.4, df09.frac_enr, yerr = df09.frac_enr_err, fmt="o", color="red", label="0.9<PIP")
plt.legend(title = "Whole Blood PIP bin", loc="upper left")#, fontsize=12)
plt.xlabel("Normalized EMS bin")
plt.ylabel("Enrichment")
plt.xticks(np.arange(6), xtk, rotation=40)
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/fig1supp_emsitself_wb.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()

