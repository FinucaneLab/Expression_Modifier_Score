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


##validation of pip_EMS with SuRE data
df1 = hl.read_table("gs://qingbowang/sure_data/k562_case_maxems_pipuf.ht").to_pandas()
df1["case"] = True
def pipbin(x):
    if x<0.01: return (-1)
    elif x<0.1: return (0)
    elif x<0.5: return (1)
    elif x<0.9: return (5)
    else: return (9)
df["u_bin"] = df.max_pip_u.apply(lambda x: pipbin(x))
df["f_bin"] = df.max_pip_f.apply(lambda x: pipbin(x))
#get the full data
fd = hl.read_table("gs://qingbowang/sure_data/fulldata_pipuf.ht")
fd = fd.filter(hl.is_defined(fd.max_pip_u) & (hl.is_defined(fd.max_pip_f)))
fdf = fd.select("k562.wilcox.p.value", "max_pip_u", "max_pip_f").to_pandas()
fdf["u_bin"] = fdf.max_pip_u.apply(lambda x: pipbin(x))
fdf["f_bin"] = fdf.max_pip_f.apply(lambda x: pipbin(x))
#annotate the SuRE positives
fdf.index = fdf.v
df1 = df1[~df1.v.isna()]
df1.index = df1.v
fdf = fdf.join(df1["case"], how="left")
fdf["case"] = fdf.case.fillna(False)

#plot
d1 = fdf.groupby(["u_bin"]).case.agg(["mean","sem", "count"])
d2 = fdf.groupby(["f_bin"]).case.agg(["mean","sem", "count"])
plt.rcParams.update({'font.size': 16})
F = sum(fdf.case)/len(fdf.case)
plt.figure(figsize=(5, 4))
plt.errorbar(np.arange(5), d1["mean"]/F, d1["sem"]/F, color="blue", fmt="o", label="PIP$_u$")
plt.errorbar(np.arange(5)+0.1, d2["mean"]/F, d2["sem"]/F, color="tab:orange", fmt="o", label="PIP$_{EMS}$")
plt.xticks(np.arange(5), tick)
plt.xlabel("PIP bin", fontsize=16)
plt.ylabel("raQTL enrichment", fontsize=16)
plt.xticks(fontsize=14, rotation=20)
plt.yticks(fontsize=14)
plt.ylim([0,20])
plt.legend()
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/fig2_raqtl_enr.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()



##validation with possible (PIP>0.1) complex trait causal variants in UKBB
ht = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/Whole_Blood_max_pip_u_and_f_and_comp_per_v.ht")
df = ht.filter((hl.is_defined(ht.max_pip_u)) & (hl.is_defined(ht.max_pip_f))).to_pandas() #hematはなかったら0とみる
df["u_bin"] = df.max_pip_u.apply(lambda x: pipbin(x))
df["f_bin"] = df.max_pip_f.apply(lambda x: pipbin(x))

df["hemat_01"] = df.max_hemat_pip>0.1
t1 = df.groupby("u_bin").hemat_01.agg(["mean", "std", "sem", "count"])
t1["n"] = t1["mean"]*t1["count"]
t2 = df.groupby("f_bin").hemat_01.agg(["mean", "std", "sem", "count"])
t2["n"] = t2["mean"]*t2["count"]

#plot enrichment:
tick = (["<0.01", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]", "0.9<"])
plt.rcParams.update({'font.size': 16})
F = sum(df.hemat_01)/len(df.hemat_01)
plt.figure(figsize=(5, 4))
plt.errorbar(np.arange(5), t1["mean"]/F, t1["sem"]/F, color="blue", fmt="o", label="PIP$_u$")
plt.errorbar(np.arange(5)+0.1, t2["mean"]/F, t2["sem"]/F, color="tab:orange", fmt="o", label="PIP$_{EMS}$")
plt.xticks(np.arange(5), tick)
plt.xlabel("PIP bin", fontsize=18)
plt.ylabel("UKBB enrichment\n(threshold: PIP>0.1)", fontsize=16)
plt.xticks(fontsize=16, rotation=20)
plt.yticks(fontsize=16)
#plt.ylim([0,20])
plt.legend()
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/fig2_ukbb_enr_uvsf.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()

#number of co-localizing variants:
t3 = pd.DataFrame({"u":t1.n, "f":t2.n})
tk = (["<0.01", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]", "0.9<"])
plt.rcParams.update({'font.size': 16})
fig, (ax1,ax2) = plt.subplots(2,1,sharex=True,
                         figsize=(6,4.6))
ax1.spines['bottom'].set_visible(False)
ax1.tick_params(axis='x',which='both',bottom=False)
ax2.spines['top'].set_visible(False)

ax2.set_ylim(0,450)
ax1.set_ylim(8150,8600) #これを等間隔にする必要..
#ax1.set_yticks(np.arange(1000,1501,100)) #if changing the y ticks

n = np.arange(5)

bars1 = ax1.bar(n, t3.u, color="blue", width=0.3, label="PIP$_u$")
bars2 = ax2.bar(n, t3.u, color="blue", width=0.3)
bars1 = ax1.bar(n+0.3, t3.f, color="tab:orange", width=0.3, label="PIP$_{EMS}$")
bars2 = ax2.bar(n+0.3, t3.f, color="tab:orange", width=0.3)

for tick in ax2.get_xticklabels():
    tick.set_rotation(90)
d = .005
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
kwargs.update(transform=ax2.transAxes)
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
fig.text(0, 0.58, "#Variants (PIP>0.1 in UKBB)", va='center', rotation='vertical', fontsize=16)
plt.xlabel("PIP bin", fontsize=16)
ax2.set_xticks(n)
ax2.set_xticklabels(tk, rotation=20)
ax1.legend(prop={"size":18})
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/fig2_comp_numbers_fandu.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()
