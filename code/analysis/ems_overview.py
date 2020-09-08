import sys
import pandas as pd
import numpy as np
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import math
from matplotlib.colors import LogNorm
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

tissue_name = "Whole_Blood"

#aggregation per TSS distance, eQTL p value, and MAC in GTEx
ems = hl.read_table("gs://qingbowang/ems_v1_test/ems_pcausal_gtexvg_all{0}.ht".format(tissue_name))
vg = hl.read_table("gs://qingbowang/ems_v1_test/{0}_allpairs.ht".format(tissue_name))
vg = vg.annotate(vg=vg.variant_id + "_" + vg.gene_id)
vg = vg.key_by("vg")

ems = ems.join(vg, how="left")
ems = ems.annotate(conf_gain_log10_bin = hl.ceil(ems.confidence_gain_log10))

#tss dist bin
ems = ems.annotate(tss_dist_bin_unsigned = hl.ceil(hl.log10(hl.abs(ems.tss_distance))))
ems = ems.transmute(tss_dist_bin = hl.cond(ems.tss_distance>0, ems.tss_dist_bin_unsigned, ems.tss_dist_bin_unsigned * -1))
agged = ems.group_by("tss_dist_bin", "conf_gain_log10_bin").aggregate(n=hl.agg.count())
agged.export("gs://qingbowang/ems_v1_test/tmp/{0}_tssdist_vs_EMS.tsv".format(tissue_name))

#p value
ems = ems.annotate(pval_bin = hl.case()
                                       .when(ems.pval_nominal<5*10**-8, -1)
                                       .when(ems.pval_nominal>0.05, 1)
                                       .default(0))
agged = ems.group_by("pval_bin", "conf_gain_log10_bin").aggregate(n=hl.agg.count())
agged.export("gs://qingbowang/ems_v1_test/tmp/{0}_pval_vs_EMS.tsv".format(tissue_name))

#minor allele count
agged = ems.group_by("ma_count", "conf_gain_log10_bin").aggregate(n=hl.agg.count())
agged.export("gs://qingbowang/ems_v1_test/tmp/{0}_mac_vs_EMS.tsv".format(tissue_name))



#read:
fn1 = "gs://qingbowang/ems_v1_test/tmp/{0}_mac_vs_EMS.tsv".format(tissue_name)
fn2 = "gs://qingbowang/ems_v1_test/tmp/{0}_pval_vs_EMS.tsv".format(tissue_name)
fn3 = "gs://qingbowang/ems_v1_test/tmp/{0}_tssdist_vs_EMS.tsv".format(tissue_name)
with hl.hadoop_open(fn1, 'r') as f:
    df1 = pd.read_csv(f, sep="\t")
with hl.hadoop_open(fn2, 'r') as f:
    df2 = pd.read_csv(f, sep="\t")
with hl.hadoop_open(fn3, 'r') as f:
    df3 = pd.read_csv(f, sep="\t")

#plot vs MAC
df1 = df1.pivot(index="ma_count",columns="conf_gain_log10_bin", values="n").fillna(0)
ytk = [1,50,100,150,200,250,300,350,400,450,500,550,600, 650]
log_norm = LogNorm(vmin=1, vmax=10**7) #10^5がmax
cbar_ticks = [math.pow(10, i) for i in range(math.floor(math.log10(1)), 1+math.ceil(math.log10(10**7)))]
plt.figure(figsize=(7,7))
sns.heatmap((df1+1), cmap="viridis", cbar_kws={'label': 'count', "shrink": .65, "ticks":cbar_ticks}, norm=log_norm)
plt.yticks(ytk,ytk, rotation=20, fontsize=16) 
plt.xticks([0.5,1.5,2.5,3.5,4.5,5.5],["<0.1","(0.1,1]","[1,10)","[10,100)","[100,1000)","1000<"], rotation=25, fontsize=16) 
plt.xlabel("Normalized EMS bin", fontsize=20)
plt.ylabel("Minor allele count", fontsize=20)
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/mac_vs_ems.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()

#plot vs marginal p-value
df2 = df2.pivot(index="pval_bin",columns="conf_gain_log10_bin", values="n").fillna(0)
ytk = ["$<5\cdot10^8$","[$5\cdot10^8$,0.05)", "0.05<"]
log_norm = LogNorm(vmin=1, vmax=10**7) #10^5がmax
cbar_ticks = [math.pow(10, i) for i in range(math.floor(math.log10(1)), 1+math.ceil(math.log10(10**7)))]
plt.figure(figsize=(7,6))
sns.heatmap((df2+1), cmap="viridis", cbar_kws={'label': 'count', "shrink": .65, "ticks":cbar_ticks}, norm=log_norm)
plt.yticks([0.5,1.5,2.5],ytk, rotation=25, fontsize=16) 
plt.xticks([0.5,1.5,2.5,3.5,4.5,5.5],["<0.1","(0.1,1]","[1,10)","[10,100)","[100,1000)","1000<"], rotation=25, fontsize=14) 
plt.xlabel("Normalized EMS bin", fontsize=18)
plt.ylabel("Association P-value bin", fontsize=18)
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/pval_vs_ems.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()

#plot vs TSS distance
df3 = df3.pivot(index="tss_dist_bin",columns="conf_gain_log10_bin", values="n").fillna(0)
df3 = df3.iloc[:-1,:]
ytk = ["[-$10^6$,-$10^5$]","(-$10^5$,-$10^4$]","(-$10^4$,-$10^3$]","(-$10^3$,-$10^2$]","(-$10^2$,-10]","(-10,-1]","1","[1,10)","[10,$10^2$)","[$10^2$,$10^3$)","[$10^3$,$10^4$)","[$10^4$,$10^5$)","[$10^5$,$10^6$]"]
log_norm = LogNorm(vmin=1, vmax=10**7) #10^5がmax
cbar_ticks = [math.pow(10, i) for i in range(math.floor(math.log10(1)), 1+math.ceil(math.log10(10**7)))]
plt.figure(figsize=(7,7))
sns.heatmap((df3+1), cmap="viridis", cbar_kws={'label': 'count', "shrink": .65, "ticks":cbar_ticks}, norm=log_norm)
plt.yticks(np.arange(13)+0.5,ytk, rotation=0, fontsize=16) 
plt.xticks([0.5,1.5,2.5,3.5,4.5,5.5],["<0.1","(0.1,1]","[1,10)","[10,100)","[100,1000)","1000<"], rotation=28, fontsize=15) 
plt.xlabel("Normalized EMS bin", fontsize=20)
plt.ylabel("Distance to TSS bin", fontsize=20)
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/tss_vs_ems.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()




