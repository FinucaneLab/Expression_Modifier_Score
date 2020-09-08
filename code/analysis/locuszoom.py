import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"


tissue_name = "Whole_Blood"

#focus on the gene to plot:
ht = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/Whole_Blood_updated_pips_alpha_filtered_by_cs.ht")
pipf = ht.filter((ht.gene=="ENSG00000109046.14")).to_pandas()

#plot
pipf["pos"] = pipf.variant_hg38.str.split("_").str[1].astype(int)
pipf["tss_dist"] = pipf.pos - 27245021 - 49055
#pipf["ems_norm"] = pipf.ems / 1.213517229204564e-05 #いや、そのままでいいや..
f, (ax1, ax2, ax3) = plt.subplots(3,1, sharex=True)

ax1.scatter(pipf[abs(pipf.tss_dist + 49055)<10000].tss_dist, pipf[abs(pipf.tss_dist + 49055)<10000].pip, color="blue", label="pip")
ax2.scatter(pipf[abs(pipf.tss_dist + 49055)<10000].tss_dist, pipf[abs(pipf.tss_dist + 49055)<10000].ems, color="green", label="EMS")
ax3.scatter(pipf[abs(pipf.tss_dist + 49055)<10000].tss_dist, pipf[abs(pipf.tss_dist + 49055)<10000].updated_pip, color="orange")
ax3.set_ylabel("PIP$_{EMS}$")
ax2.set_ylabel("EMS")
ax1.set_ylabel("PIP$_u$")
ax3.set_xlabel("Distance to TSS")
ax1.set_ylim([0,1.05])
ax3.set_ylim([0,1.05])
ax1.axvline(x=-49055, linestyle="--", linewidth=0.5, color="black")
ax2.axvline(x=-49055, linestyle="--", linewidth=0.5, color="black")
ax3.axvline(x=-49055, linestyle="--", linewidth=0.5, color="black")

plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/fig4_example_plot_01.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()

#zoomed out plot:
pipf.sort_values(by="pip", inplace=True, ascending=False)
pipf["cum_pip"] = np.cumsum(pipf.pip)
thres = pipf[pipf.cum_pip>0.95].index[0]
pipf_cs = pipf.iloc[:thres+1,:]
pipf_else = pipf.iloc[thres+1:,:]

f, (ax1, ax2, ax3) = plt.subplots(3,1, sharex=True)

ax1.scatter(pipf_cs.tss_dist, pipf_cs.pip, color="blue")
ax2.scatter(pipf_cs.tss_dist, pipf_cs.ems, color="green")
ax3.scatter(pipf_cs.tss_dist, pipf_cs.updated_pip, color="orange")
ax1.scatter(pipf_else.tss_dist, pipf_else.pip, color="grey", label="outside the 95% CS")
ax2.scatter(pipf_else.tss_dist, pipf_else.ems, color="grey")
ax3.scatter(pipf_else.tss_dist, pipf_else.updated_pip, color="grey")
ax1.legend()

ax3.set_ylabel("PIP$_{EMS}$")
ax2.set_ylabel("EMS")
ax1.set_ylabel("PIP$_u$")
ax3.set_xlabel("Distance to TSS")
ax1.axvline(x=-49055, linestyle="--", linewidth=0.5, color="black")
ax2.axvline(x=-49055, linestyle="--", linewidth=0.5, color="black")
ax3.axvline(x=-49055, linestyle="--", linewidth=0.5, color="black")

plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/fig4_example_plot_macro_color.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()
