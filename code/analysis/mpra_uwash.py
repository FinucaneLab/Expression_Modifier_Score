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


#load the results:
fn = "gs://qingbowang/ems_v1_test/fig2/ems_mpra_table.tsv"
with hl.hadoop_open(fn, 'r') as f:
    df = pd.read_csv(f, sep="\t", index_col=0)
fn = "gs://qingbowang/tmp/uwash_mpra_ncer.tsv"
with hl.hadoop_open(fn, 'r') as f:
    nc = pd.read_csv(f, sep="\t")
nc.head()
#concat the results:
df.index = df.iloc[:,0] + "_" + df.ID
nc.index = nc["#CHROM"]+"_"+nc.POS.astype(str)+"_"+nc.REF+"_"+nc.ALT+"_"+nc.ID
df["ncer"] = nc.loc[df.index, "ncer"]

#plot
df.bonf_signif = (df.pval<0.05/df.shape[0])
denom = df.bonf_signinf.sum() / df.shape[0]
plt.rcParams.update({'font.size': 16})
fig, ax0 = plt.subplots(1, figsize=(10,4))
colors = ["tab:green", "tab:blue", "tab:pink", "tab:olive", "tab:cyan", "tab:brown"]

i = 0
for s in ["rawp_gmean", "deepsea_e","cadd","ncer", "gerp","fathmm"]:
    df["x_bin"] = df[s].rank(pct=True) // 0.1
    y_bin = df.groupby("x_bin").bonf_signinf.agg(['mean','sem'])
    ax0.errorbar((np.arange(10)-5)/20 + i, y_bin["mean"]/denom, yerr=y_bin["sem"]/denom, fmt="o", color=colors[i])
    if i==0:
        ax0.axhline(y=y_bin["mean"].values[0]/denom, color="green", linestyle="--", linewidth=.5)
        ax0.axhline(y=y_bin["mean"].values[-1]/denom, color="green", linestyle="--", linewidth=.5)    
    i += 1
    #print (s)
    #print (y_bin/denom)
#ax0.set_title("enhancer variants, n={0}".format(df.shape[0]))
plt.xticks(np.arange(6), ["EMS","DeepSEA","CADD", "ncER", "GERP","Fathmm"], size=20, rotation=20)
#ax0.set_title("promoter + enhancer variants, n={0}".format(df.shape[0]))
plt.xlabel("Score (decile)", size=20)
plt.ylabel("MPRA hits enrichment", size=18)
fig.tight_layout()
plt.subplots_adjust(hspace=0.4)
fn = "gs://qingbowang/ems_v1_test/fig2/fig2_mpra_combined_enr_wncer.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()

#and per element:
tab = {}
for s in ["rawp_gmean", "deepsea_e","fathmm","cadd","gerp", "ncer"]:
    tab[s] = []
    for ele in df.target_gene.unique():
        df1 = df[df.target_gene==ele]
        df1["x_bin"] = abs(df1[s]).rank(ascending=False)
        y_m = df1[df1["x_bin"]<df1.shape[0]/10].bonf_signinf.sum() #top 10%
        y_all = df1.bonf_signinf.sum()
        tab[s].append((y_m/y_all))
tab["n"] = []
tab["n_pos"] = []
tab["enh"] = []
for ele in df.target_gene.unique():
    tab["n"].append(sum(df.target_gene==ele))
    tab["n_pos"].append(sum((df.target_gene==ele) & (df.bonf_signinf)))
    tab["enh"].append(df[df.target_gene==ele].enh.values[0])
tab = pd.DataFrame(tab)
tab.index = df.target_gene.unique()
tab = tab[tab.n_pos>50]
tab.sort_values(by=["enh", "n_pos"])
enh = tab[tab.enh=="enh"].sort_index().sort_values(by="enh").iloc[:,:-3]
prom = tab[tab.enh=="prom"].sort_index().sort_values(by="enh").iloc[:,:-3]
w = pd.DataFrame([np.nan]*6)
w.index = enh.columns
w.columns = [""]
w = w.T
toplot = pd.concat([enh, w, prom])

maxlabel = toplot.apply(lambda x: x.rank(), axis=1)
maxlabel[maxlabel!=5] = ""
maxlabel[maxlabel==5] = "*"
#match the order
toplot = toplot.iloc[:,[0,1,3,5,4,2]]
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"


toplot.columns = ["EMS","DeepSEA","CADD", "ncER", "GERP","Fathmm"]
plt.figure(figsize = (5,8))
sns.heatmap(toplot*10**2, cmap="cool", xticklabels=True, yticklabels=True, annot = maxlabel, fmt="")
plt.title("% positive hits in top decile\n *=best score")
plt.xlabel("score")
plt.ylabel("functional element\n (promoter)                    (enhancer)")
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/fig2/fig2_mpra_supp_elements_w_ncer.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()



#mean vs median etc
f = "gs://qingbowang/ems_validation_data/uwash_mpra_ems_confidence_gain_rest20tissues_fixedtssdist.tsv"
ems1 = hl.import_table(f, impute=True).to_pandas() #This is the baseline score
ems1.index = ems1[""]
del ems1[""]
f = "gs://qingbowang/ems_validation_data/uwash_mpra_ems_confidence_gain_30tissues_fixedtssdist.tsv"
ems2 = hl.import_table(f, impute=True).to_pandas() #This is the baseline score
ems2.index = ems2[""]
del ems2[""]
ems = pd.concat([ems1, ems2], axis=1)
ems.index = ems.index.str.replace("uwash_mpra_target_","")


df.index = df.index + "_" + df.target_gene

mat = pd.concat([t[["rawp_mean","rawp_median","rawp_gmean","rawp_min","rawp_max"]], t.iloc[:,-49:]], axis=1) #matrix to plot
mat.columns = mat.columns.str.replace("rawp_","")
label = pd.DataFrame(t.bonf_signinf)
ys = {}
yerrs = {}
def decile_enr(v, label):
    rk = v.rank(pct=True)//0.1
    label["x_bin"] = rk
    y_bin = label.groupby("x_bin").bonf_signinf.agg(['mean','sem'])
    return (y_bin["mean"],y_bin["sem"])
for s in mat.columns:
    (y,yerr) = decile_enr(mat[s], label)
    ys[s] = y
    yerrs[s] = yerr
ys = pd.DataFrame(ys)
yerrs = pd.DataFrame(yerrs)
#plot
fig, (ax0) = plt.subplots(1, sharey=True, sharex=True, figsize=(12,3))
i = 0
for i in range(5):
    ax0.errorbar((np.arange(10)-5)/20 + i, ys.iloc[:,i], yerr=yerrs.iloc[:,i], fmt="o", color="green")
    if i==2:
        ax0.axhline(y=ys.iloc[:,i].values[0], color="green", linestyle="--", linewidth=.5)
        ax0.axhline(y=ys.iloc[:,i].values[-1], color="green", linestyle="--", linewidth=.5)    

    i += 1
ax0.set_ylim([0,0.45])
plt.xticks(np.arange(5), ys.columns[:5], size=18)
plt.xlabel("method to combine EMS (decile)", size=18)
plt.ylabel("fraction of MPRA hits", size=16)
fig.tight_layout()
plt.subplots_adjust(hspace=0.4)
fn = "gs://qingbowang/ems_v1_test/fig2/figs_mpra_meanvsmedian.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()


#per tissue:
fn = "gs://qingbowang/gtex_colors.tsv"
with hl.hadoop_open(fn, 'r') as f:
        gtex_colors = pd.read_csv(f, sep="\t", index_col=0)
gtex_colors["tissue_color_hex"] = "#" + gtex_colors["tissue_color_hex"].astype(str) 
fn = "gs://qingbowang/gtex_name_corres.csv"
with hl.hadoop_open(fn, 'r') as f:
        gtex_names = pd.read_csv(f, sep=",", index_col=0)
gtex_colors = pd.concat([gtex_colors, gtex_names], axis=1).fillna("0")
gtex_colors.index = gtex_colors.my_name.str.replace(" ", "")
ys.columns = ys.columns.str.replace("\n","_") #putting this back

fig, axes = plt.subplots(5,1, sharey=True, figsize=(15,12)) #sharexしてるからaxisうまくいかんぽい.
j = 0 #j行目
i = 0
for i in range(10):
    cl = gtex_colors.loc[ys.columns[5+i+j*10],"tissue_color_hex"]#color
    axes[0].errorbar((np.arange(10)-5)/15 + i, ys.iloc[:,5+i+j*10], yerr=yerrs.iloc[:,5+i+j*10], fmt="o", color=cl)
    i += 1
axes[0].set_ylim([0,0.5])
plt.sca(axes[0])
plt.xticks(np.arange(10), gtex_colors.loc[ys.columns[j*10+5:(j+1)*10+5],"tissue_abbrv"], rotation=0, size=14)

j = 1
for i in range(10):
    cl = gtex_colors.loc[ys.columns[5+i+j*10],"tissue_color_hex"]#color
    axes[1].errorbar((np.arange(10)-5)/15 + i, ys.iloc[:,5+i+j*10], yerr=yerrs.iloc[:,5+i+j*10], fmt="o", color=cl)
    i += 1
ax1.set_ylim([0,0.5])
axes[1].set_xticks(np.arange(10))
axes[1].set_xticklabels(gtex_colors.loc[ys.columns[j*10+5:(j+1)*10+5],"tissue_abbrv"], rotation=0, size=14)

j = 2
for i in range(10):
    cl = gtex_colors.loc[ys.columns[5+i+j*10],"tissue_color_hex"]#color
    axes[2].errorbar((np.arange(10)-5)/15 + i, ys.iloc[:,5+i+j*10], yerr=yerrs.iloc[:,5+i+j*10], fmt="o", color=cl)
    i += 1
axes[2].set_ylim([0,0.5])
axes[2].set_xticks(np.arange(10))
axes[2].set_xticklabels(gtex_colors.loc[ys.columns[j*10+5:(j+1)*10+5],"tissue_abbrv"], rotation=0, size=14)

j = 3
for i in range(10):
    cl = gtex_colors.loc[ys.columns[5+i+j*10],"tissue_color_hex"]#color
    axes[3].errorbar((np.arange(10)-5)/15 + i, ys.iloc[:,5+i+j*10], yerr=yerrs.iloc[:,5+i+j*10], fmt="o", color=cl)
    i += 1
axes[3].set_ylim([0,0.5])
axes[3].set_xticks(np.arange(10))
axes[3].set_xticklabels(gtex_colors.loc[ys.columns[j*10+5:(j+1)*10+5],"tissue_abbrv"], rotation=0, size=14)

j = 4
for i in range(9):
    cl = gtex_colors.loc[ys.columns[5+i+j*10],"tissue_color_hex"]#color
    axes[4].errorbar((np.arange(10)-5)/15 + i, ys.iloc[:,5+i+j*10], yerr=yerrs.iloc[:,5+i+j*10], fmt="o", color=cl)
    i += 1
#dummy for i=10
axes[4].errorbar((np.arange(10)-5)/15 + i, ys.iloc[:,5+i-1+j*10], yerr=yerrs.iloc[:,5+i-1+j*10], fmt="o", color="#00000000")
axes[4].set_ylim([0,0.5])
axes[4].set_xticks(np.arange(10))
axes[4].set_xticklabels(list(gtex_colors.loc[ys.columns[j*10+5:(j+1)*10+5],"tissue_abbrv"]) + [""], rotation=0, size=14)

axes[4].set_xlabel("tissue to train EMS (decile)", size=18)
axes[2].set_ylabel("fraction of MPRA hits", size=18)
fig.tight_layout()
plt.subplots_adjust(hspace=0.4)
fn = "gs://qingbowang/ems_v1_test/fig2/figs_mpra_pertissue.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=350)
plt.show()


