import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')
import pandas as pd
import numpy as np

tissue_name = "Whole_Blood"
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

#aggregation:
ht = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_alpha_filtered_by_cs.ht".format(tissue_name))
ht = ht.annotate(pip_bin_unif_prior = hl.case()
                                             .when(ht.pip > 0.9, 9)
                                             .when(ht.pip > 0.5, 5)
                                             .when(ht.pip > 0.1, 1)
                                             .when(ht.pip > 0.01, 0)
                                             .default(-1),
                 pip_bin_ems_as_a_prior = hl.case()
                 .when(ht.updated_pip > 0.9, 9)
                 .when(ht.updated_pip > 0.5, 5)
                 .when(ht.updated_pip > 0.1, 1)
                 .when(ht.updated_pip > 0.01, 0)
                 .default(-1)
                 )
ht.group_by("pip_bin_unif_prior", "pip_bin_ems_as_a_prior").aggregate(n=hl.agg.count()).export("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_overview_cutoff_cs.tsv".format(tissue_name))
ht.filter(ht.updated_pip>0.00425).write("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_csfiltered_00425.ht".format(tissue_name)) #for more inspection
ht.group_by("gene","cs_id").aggregate(n=hl.agg.count(), M_u = hl.agg.max(ht.pip), M_f=hl.agg.max(ht.updated_pip)).write("gs://qingbowang/ems_v1_test/updated_pips/{0}_pip_u_cs_basic_stats_csfiltered.ht".format(tissue_name))



#count the number of variants per CS
df = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/Whole_Blood_pip_u_cs_basic_stats_csfiltered.ht").to_pandas()
dff = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/Whole_Blood_updated_pips_csfiltered_00425.ht").to_pandas()

df["n_u"] = df.n

def get_n_f(g, cs_id):
    v = dff[(dff.gene==g) & (dff.cs_id==cs_id)].sort_values(by="updated_pip", ascending=False)
    if v.shape[0]==0:
        return (1000000) #just a sign for "more than 20"
    else:
        v["cum"] = v.updated_pip.cumsum()
        n = v.shape[0] - sum(v.cum>0.95) + 1
        return (n)
df["n_f"] = df.apply(lambda x: get_n_f(x["gene"] ,x["cs_id"]), axis=1)

#bin it:
def bin_cnt(x):
    if x<5:
        return str(x)
    elif x<11:
        return ("[5,10]")
    elif x<21:
        return ("[11,20]")
    else:
        return ("20<")
df["n_u_bin"] = df.n_u.apply(lambda x: bin_cnt(x))
df["n_f_bin"] = df.n_f.apply(lambda x: bin_cnt(x))
table = df[["n_u_bin", "n_f_bin"]]
#save;
fn = "gs://qingbowang/ems_v1_test/updated_pips/Whole_Blood_n_var_per_cs_comparison_csfiltered.tsv"
hl.Table.from_pandas(table).export(fn)
#read 
fn = "gs://qingbowang/ems_v1_test/updated_pips/Whole_Blood_n_var_per_cs_comparison_csfiltered.tsv"
df = hl.import_table(fn).to_pandas()

#make it a table for visualization:
table = df.groupby(["n_u_bin", "n_f_bin"]).size().unstack().fillna(0).astype(int)
table = table[["1","2","3","4","[5,10]","[11,20]","20<"]]
table = table.T
table = table[["1","2","3","4","[5,10]","[11,20]","20<"]]
table = table.T

#plot
import math
from matplotlib.colors import LogNorm
log_norm = LogNorm()#vmin=table.min().min()+1, vmax=table.max().max()+1)
plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(10,10))
sns.heatmap(table+1, annot=table, fmt="d", square=True, linewidths=.5, norm=log_norm,
            cmap="viridis", cbar_kws={'label': 'count', 
                                      "shrink": .6})#,
                                      #"ticks":{0:0,1:2,2:4,3:8}})#, vmin=0, vmax=3000)
plt.yticks(rotation=20, fontsize=20) 
plt.xticks(rotation=20, fontsize=20) 
plt.xlabel("#variants in CS using EMS as a prior", fontsize=22)
plt.ylabel("#variants in CS using uniform prior", fontsize=22)
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/fig4_comp_n_u_vs_f_csfiltered_v2.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()

#plot the per variant table
fn = "gs://qingbowang/ems_v1_test/updated_pips/Whole_Blood_updated_pips_overview_cutoff_cs.tsv"
with hl.hadoop_open(fn, 'r') as f:
    df = pd.read_csv(f, sep="\t")

df = df.pivot(index='pip_bin_unif_prior', columns='pip_bin_ems_as_a_prior', values='n').fillna(0).astype(int)
idx = ["<0.01","(0.01,0.1]","(0.1,0.5]","(0.5,0.9]","0.9<"]
df.index = idx
df.columns = idx
df.index.names = ["pip with uniform prior"]
df.columns.names = ["updated pip using EMS"]
df_rownorm = df.apply(lambda x: x / sum(x), axis=1)
log_norm = LogNorm(vmin=1, vmax=10**5) #10^5がmax
cbar_ticks = [math.pow(10, i) for i in range(math.floor(math.log10(1)), 1+math.ceil(math.log10(10**5)))]

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(10,10))
sns.heatmap((df+1).applymap(lambda x: min(10**5, x)), annot=df, fmt="d", square=True, linewidths=.5, norm=log_norm, 
            cmap="viridis", cbar_kws={'label': 'count', "shrink": .65, "ticks":cbar_ticks})#, vmin=0, vmax=100)
plt.yticks(rotation=20, fontsize=20) 
plt.xticks(rotation=20, fontsize=20) 
plt.xlabel("PIP using EMS as a prior", fontsize=24)
plt.ylabel("PIP using uniform prior", fontsize=24)
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/fig4_pip_u_vs_f.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()
