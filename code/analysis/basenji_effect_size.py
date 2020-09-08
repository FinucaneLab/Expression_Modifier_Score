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

#the basenji cage score:
pos = hl.import_table("gs://qingbowang/ems_v1_test/basenji_cage_neutrophil_wb_pos.tsv", impute=True)
#the effect size in gtex:
posd = hl.import_table("gs://qingbowang/ems_v1_test/{0}_positive_training_vg_beforeannot.tsv".format(tissue_name), impute=True)
posd = posd.select("variant_id", "gene_id", "slope")
posd = posd.to_pandas()
posbas = hl.import_table("gs://qingbowang/ems_v1_test/basenji_cage_neutrophil_wb_pos.tsv", impute=True).to_pandas()
posd.index = posd.variant_id.str.replace("_b38","")
posbas.index = posbas.ID
df = posd.join(posbas, how="left")


#bin the scores:
def bas_bin(x):
    pos = (x>0)
    ab = abs(x)
    if ab>10: b = 10
    elif ab>1: b = 1
    elif ab>0.1: b=0.1
    else: b = 0.01
    if pos: return (b)
    else: return (b*-1)
df["bas_bin"] = df.CNhs10862.apply(lambda x: bas_bin(x))

def beta_bin(x):
    pos = (x>0)
    ab = abs(x)
    if ab>0.5: b = 5
    else: b = 1
    if pos: return (b)
    else: return (b*-1)
df["beta_bin"] = df.slope.apply(lambda x: beta_bin(x))

#get the summary:
tb = df.groupby(["bas_bin","beta_bin"]).size()
tb = pd.DataFrame(tb).unstack()
p0 = tb.sum(axis=1)
p0 = p0 / p0.sum()
p1 = tb/tb.sum(axis=0)
enr = (p1.T / np.array(p0)).T.fillna(0)
err = np.sqrt(p1*(1-p1)/tb.sum(axis=0))
err = (err.T / np.array(p0)).T.fillna(0)

#plot
plt.rcParams.update({'font.size': 20})
plt.figure(figsize=(10, 4))
plt.errorbar(np.arange(10), enr.iloc[:,0], err.iloc[:,0],label="β<-0.8", fmt="o", color="mediumblue")
plt.errorbar(np.arange(10)+15, enr.iloc[:,1], err.iloc[:,1], label="-0.8<β<0", fmt="o", color="darkturquoise")
plt.errorbar(np.arange(10)+15*2, enr.iloc[:,2], err.iloc[:,2], label="0<β<0.8", fmt="o", color="tab:pink")
plt.errorbar(np.arange(10)+15*3, enr.iloc[:,3], err.iloc[:,3], label="0.8<β", fmt="o", color="tab:red")
plt.legend(bbox_to_anchor=[1.02, 0.8])
#plt.xticks([5,20,35,50])
plt.xticks([])
plt.xlabel("Basenji score bin for each β")
plt.ylabel("Enrichment")
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/bas_eff_size.png"
with hl.hadoop_open(fn, 'wb') as f:
            plt.savefig(f, dpi=400)

plt.show()

