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

vg = hl.read_table("gs://qingbowang/ems_v1_test/{0}_allpairs_binannot_fmannot.ht".format(tissue_name))
vg.head(5).show()

vg = vg.annotate(susie_bin = hl.case()
                                      .when(vg.pp_susie>0.9, 9)
                                      .when(vg.pp_susie > 0.5, 5)
                                      .when(vg.pp_susie > 0.1, 1)
                                      .when(vg.pp_susie > 0.01, 0)
                                      .default(-1),
                 fm_bin = hl.case()
                                      .when(vg.pp_fm>0.9, 9)
                                      .when(vg.pp_fm > 0.5, 5)
                                      .when(vg.pp_fm > 0.1, 1)
                                      .when(vg.pp_fm > 0.01, 0)
                 .default(-1))
#aggregation:
for f in ["Intron_UCSC", "Enhancer_Hoffman","DHS_Trynka","TFBS_ENCODE", "UTR_3_UCSC", "Promoter_UCSC", "Coding_UCSC","UTR_5_UCSC"]:
    vg.group_by(f, "susie_bin","fm_bin").aggregate(n = hl.agg.count()).export("gs://qingbowang/ems_v1_test/{0}_{1}_fm_and_susie_bin.tsv".format(tissue_name, f))
    
#read as tsv
enrs = {}
errs = {}
for feat in ["Intron_UCSC", "Enhancer_Hoffman","DHS_Trynka","TFBS_ENCODE", "UTR_3_UCSC", "Promoter_UCSC", "Coding_UCSC","UTR_5_UCSC"]:
    fn = "gs://qingbowang/ems_v1_test/{0}_{1}_fm_and_susie_bin.tsv".format(tissue_name, feat)
    with hl.hadoop_open(fn, 'r') as f:
            df = pd.read_csv(f, sep="\t")
    dfpos = df[df[feat]].pivot(index='susie_bin', columns='fm_bin', values='n')#n(positive)
    dfneg = df[~df[feat]].pivot(index='susie_bin', columns='fm_bin', values='n')#n(negative)
    enr = dfpos/(dfpos+dfneg)
    err = np.sqrt(enr*(1-enr)/(dfpos+dfneg))
    denom = dfpos.sum().sum()/(dfpos.sum().sum()+dfneg.sum().sum())
    enr = enr/denom
    err = err/denom
    enrs[feat] = enr
    errs[feat] = err

#plot enhancer:
plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(9, 4.5))
tick = (["<0.01", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]", "0.9<"])
toplot = enrs["Enhancer_Hoffman"]
toplot_err = errs["Enhancer_Hoffman"]
plt.errorbar(np.arange(5), toplot.iloc[0,:], toplot_err.iloc[0,:], color="grey", fmt="o", label="FINEMAP PIP<0.01")
plt.errorbar(np.arange(5)+0.1, toplot.iloc[1,:], toplot_err.iloc[1,:], color="khaki", fmt="o", label="0.01<FINEMAP PIP<0.1")
plt.errorbar(np.arange(5)+0.2, toplot.iloc[2,:], toplot_err.iloc[2,:], color="gold", fmt="o", label="0.1<FINEMAP PIP<0.5")
plt.errorbar(np.arange(5)+0.3, toplot.iloc[3,:], toplot_err.iloc[3,:], color="orange", fmt="o", label="0.5<FINEMAP PIP<0.9")
plt.errorbar(np.arange(5)+0.4, toplot.iloc[4,:], toplot_err.iloc[4,:], color="red", fmt="o", label="FINEMAP PIP 0.9<")

plt.legend(bbox_to_anchor= (1.01, 1.01))
plt.xlim([-0.1,4.5])
plt.xticks(np.arange(5), tick, rotation=30)
plt.xlabel("SuSiE PIP bin", size=16)
plt.ylabel("Enhancer enrichment", size=18)
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/susie_vs_fm_enh.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()

#plot intron:
plt.figure(figsize=(5.5, 4.5))
tick = (["<0.01", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]", "0.9<"])
toplot = enrs["Intron_UCSC"]
toplot_err = errs["Intron_UCSC"]
plt.errorbar(np.arange(5), toplot.iloc[0,:], toplot_err.iloc[0,:], color="grey", fmt="o", label="FINEMAP PIP<0.01")
plt.errorbar(np.arange(5)+0.1, toplot.iloc[1,:], toplot_err.iloc[1,:], color="khaki", fmt="o", label="0.01<FINEMAP PIP<0.1")
plt.errorbar(np.arange(5)+0.2, toplot.iloc[2,:], toplot_err.iloc[2,:], color="gold", fmt="o", label="0.1<FINEMAP PIP<0.5")
plt.errorbar(np.arange(5)+0.3, toplot.iloc[3,:], toplot_err.iloc[3,:], color="orange", fmt="o", label="0.5<FINEMAP PIP<0.9")
plt.errorbar(np.arange(5)+0.4, toplot.iloc[4,:], toplot_err.iloc[4,:], color="red", fmt="o", label="0.9<FINEMAP PIP")

plt.xlim([-0.1,4.5])
plt.xticks(np.arange(5), tick, rotation=30)
plt.xlabel("SuSiE PIP bin", size=16)
plt.ylabel("Intron enrichment", size=18)
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/susie_vs_fm_intron.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()
