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


tissues = ["Whole_Blood",
            "Muscle_Skeletal",
           "Liver",
           "Brain_Cerebellum",
            "Prostate",
            "Spleen",
            "Skin_Sun_Exposed_Lower_leg",
            "Artery_Coronary",
            "Esophagus_Muscularis",
            "Esophagus_Gastroesophageal_Junction",
            "Artery_Tibial",
            "Heart_Atrial_Appendage",
            "Nerve_Tibial",
            "Heart_Left_Ventricle",
            "Adrenal_Gland",
            "Adipose_Visceral_Omentum",
            "Pancreas",
            "Lung",
            "Pituitary",
            "Brain_Nucleus_accumbens_basal_ganglia",
            "Colon_Transverse",
            "Adipose_Subcutaneous",
            "Esophagus_Mucosa",
            "Brain_Cortex",
            "Thyroid",
            "Stomach",
            "Breast_Mammary_Tissue",
            "Colon_Sigmoid",
            "Skin_Not_Sun_Exposed_Suprapubic",
            "Testis",
            "Artery_Aorta",
            "Brain_Amygdala",
            "Brain_Anterior_cingulate_cortex_BA24",
            "Brain_Caudate_basal_ganglia",
            "Brain_Cerebellar_Hemisphere",
            "Brain_Frontal_Cortex_BA9",
            "Brain_Hippocampus",
            "Brain_Hypothalamus",
            "Brain_Putamen_basal_ganglia",
            "Brain_Spinal_cord_cervical_c-1",
            "Brain_Substantia_nigra",
            "Cells_Cultured_fibroblasts",
            "Cells_EBV-transformed_lymphocytes",
            "Kidney_Cortex",
            "Minor_Salivary_Gland",
            "Ovary",
            "Small_Intestine_Terminal_Ileum",
            "Uterus",
            "Vagina"]

#for the color
fn = "gs://qingbowang/gtex_colors.tsv"
with hl.hadoop_open(fn, 'r') as f:
        gtex_colors = pd.read_csv(f, sep="\t", index_col=0)
gtex_colors["tissue_color_hex"] = "#" + gtex_colors["tissue_color_hex"].astype(str) 
fn = "gs://qingbowang/gtex_name_corres.csv"
with hl.hadoop_open(fn, 'r') as f:
        gtex_names = pd.read_csv(f, sep=",", index_col=0)
gtex_colors = pd.concat([gtex_colors, gtex_names], axis=1).fillna("0")
gtex_colors.index = gtex_colors.my_name.str.replace(" ", "")
            
#number of putative causal eQTLs per tissue:
npos = []
for tissue_name in tissues:
    fn = "gs://qingbowang/ems_v1_test/{0}_positive_training_vg_beforeannot.tsv".format(tissue_name)
    with hl.hadoop_open(fn, 'r') as f:
            df = pd.read_csv(f, sep="\t")
    npos.append(df.shape[0])
df = pd.DataFrame({"tissue":tissues, "n":npos})

df.index = df.tissue
df = pd.concat([df,gtex_colors.loc[df.index, :]], axis=1)

#plot:
df.sort_values(by="n", ascending=False, inplace=True)
sns.set()
sns.set(rc={'figure.figsize':(10,5), "xtick.bottom":True, "ytick.left":True, 
            'axes.facecolor': 'white', "axes.edgecolor": "black", 
            'grid.color': '#aaaaaaff', "font.sans-serif":"Arial"})
sns.barplot(x="tissue", y="n", data=df, palette=df.tissue_color_hex)
plt.xticks(np.arange(df.shape[0]), df.tissue_abbrv, rotation=90)
plt.xlabel("Tissue", size=16)
plt.ylabel("Number of putative \n causal variant-gene pairs", size=16)
plt.xticks(size=14)
plt.yticks(size=16)
fn = "gs://qingbowang/ems_v1_test/fig/fig1_a_n_causal.png"
plt.tight_layout()
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=500)
plt.show()


#the n(unique causal variant), n unique and n shared 
vgs = {}
for tissue_name in tissues:
    fn = "gs://qingbowang/ems_v1_test/{0}_positive_training_vg_beforeannot.tsv".format(tissue_name)
    with hl.hadoop_open(fn, 'r') as f:
            dfs = pd.read_csv(f, sep="\t")
    vg = dfs.variant_id + "_" + dfs.gene_id
    vgs[tissue_name] = list(vg)
vs = []
for tissue_name in tissues:
    v = pd.Series(vgs[tissue_name]).str.split("_b38").str[0]
    vs = np.union1d(vs, v)
vglong = []
for tissue_name in tissues:
    vglong = vglong + (vgs[tissue_name])
vc = pd.Series(vglong).value_counts()
vc = vc.value_counts()
vc = vc.sort_index()
vceasy = vc[vc.index<10] #bin things
vceasy.loc[10] = vc[(vc.index>9) & (vc.index<20)].sum()
vceasy.loc[20] = vc[(vc.index>19) & (vc.index<30)].sum()
vceasy.loc[30] = vc[(vc.index>29) & (vc.index<40)].sum()
vceasy.loc[40] = vc[(vc.index>39)].sum()

#plot the number of shared tissues:
plt.rcParams.update({'font.size': 20})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
vceasy.plot.bar(figsize=(8,4))
plt.grid(b=None, axis='x')
plt.xlabel("Number of shared tissues", fontsize=16)
plt.ylabel("Number of putative \ncausal variant-gene pairs", fontsize=14)
plt.xticks(np.arange(len(vceasy)), [1,2,3,4,5,6,7,8,9,"[10,20)", "[20,30)", "[30,40)","[40,50)"], rotation=90, fontsize=16)
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/fig1_supp_n_shared.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=500)
plt.show()

#cumulative:
vccum = vc.cumsum()
vccum = vccum/vc.sum()
vccum[0] = 0
vccum.sort_index(inplace=True)
vccum.plot(marker="o")
plt.xlim([0,50])
plt.ylim([0,1])
plt.grid(b=None, axis='x')
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlabel("Number of shared tissues", family="Arial", fontsize=20)
plt.ylabel("Numulative count of putative \ncausal variant-gene pairs", fontsize=18)
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/fig1_supp_n_shared_cum.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=500)
plt.show()

#number and fraction of unique per tissue:
vglong = []
for tissue_name in tissues:
    vglong = vglong + (vgs[tissue_name])
vc = pd.Series(vglong).value_counts()
nonuni = vc[vc>1].index
uni = {}
for tissue_name in tissues:
    uni[tissue_name] = len(np.setdiff1d(vgs[tissue_name], nonuni))
percuni = {}
for tissue_name in tissues:
    percuni[tissue_name] = uni[tissue_name]/len(vgs[tissue_name])
percuni = pd.Series(percuni).sort_index()
percuni = pd.concat([percuni,gtex_colors.loc[df.index, :]], axis=1)
percuni.sort_values(by=0, inplace=True, ascending=False)
percuni["x"] = np.arange(percuni.shape[0])
plt.figure(figsize=(10, 5))
plt.scatter(percuni.x, percuni[0], color=list(percuni.tissue_color_hex))
plt.grid('x', linestyle='--')
plt.grid(b=None, axis='y')
plt.xticks(percuni.x, list(percuni.tissue_abbrv), rotation=90)
plt.yticks(fontsize=14)
plt.xlabel("Tissue", fontsize=16)
plt.ylabel("Fraction of tissue-specific  \nputative causal variant-gene pairs", fontsize=14)
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/fig1_supp_n_unique.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=500)
plt.show()


## Tissue specificity of EMS

#table for causality
posall = {}
for tissue_name in tissues:
    fn = "gs://qingbowang/ems_v1_test/{0}_positive_training_vg_beforeannot.tsv".format(tissue_name)
    with hl.hadoop_open(fn, 'r') as f: 
        df = pd.read_csv(f, sep="\t")
    posall[tissue_name] = np.array(df.variant_id + "_" + df.gene_id)
    print ("done {0}".format(tissue_name))
uniall = np.array([])
for tissue_name in tissues:
    uniall = np.union1d(uniall, posall[tissue_name])
tab1 = pd.DataFrame(0, index=uniall, columns=tissues)
for tissue_name in tissues:
    tab1.loc[posall[tissue_name],tissue_name] = 1
#export:
fn = "gs://qingbowang/ems_v1_test/causal_vgs_table.tsv"
with hl.hadoop_open(fn, 'w') as f: 
    tab1.to_csv(f, sep="\t")
#table for EMS
for tissue_name in tissues:
    ht = hl.read_table("gs://qingbowang/ems_v1_test/ems_pcausal_gtexvg_all{0}.ht".format(tissue_name))
    ht = ht.key_by("vg").select("confidence_gain")
    htsub = ht[tab1ht.key]
    htsub.export("gs://qingbowang/ems_v1_test/tmp/conf_gain_for_causalvg_{0}.tsv".format(tissue_name)) #and pick it up later
conf_gain = []
for tissue_name in tissues:
    fn = "gs://qingbowang/ems_v1_test/tmp/conf_gain_for_causalvg_{0}.tsv".format(tissue_name) #actually tsv but is ht by mistake
    with hl.hadoop_open(fn, 'r') as f:
        piece = pd.read_csv(f, sep="\t", index_col=0)
    piece.columns = [tissue_name]
    piece[tissue_name] = piece[tissue_name].str.split(":").str[1].str[:-1].astype(float) #dealing with those random {:} stuff
    conf_gain.append(piece)
    print ("done {0}".format(tissue_name))
conf_gain = pd.concat(conf_gain, axis=1)

#get the tissue enrichment 
enr_med = pd.DataFrame(0, index=tissues, columns=tissues)
for i in tissues:
    for j in tissues:
        causal_indicator = tab1.loc[:,i]
        x1 = conf_gain.loc[causal_indicator==1,j]
        x2 = conf_gain.loc[causal_indicator==0,j]
        enr_med.loc[i,j] = np.log2(x1.median()/x2.median())
    print ("done i={0}".format(i))

enr_med.sort_index(inplace=True)
enr_med = enr_med.T.sort_index().T
colors = list(gtex_colors.loc[enr_med.index,"tissue_color_hex"])
enr_med.index = list(gtex_colors.loc[enr_med.index,"tissue_abbrv"])
enr_med.columns = list(enr_med.index)

#focus on whole blood:
l = pd.DataFrame(enr_med.loc["WHLBLD",:].sort_values(ascending=False))
gtex_colors.index = gtex_colors.tissue_abbrv
colors = list(gtex_colors.loc[l.index,"tissue_color_hex"])
#plot the enrichment of EMS_Whole_Blood
l["x"] = np.arange(l.shape[0])
plt.figure(figsize=(10, 5))
plt.scatter(l.x, l.WHLBLD, color = colors)
#plt.grid('x', linestyle='--')
#plt.grid(b=None, axis='y')
plt.xticks(l.x, list(l.index), rotation=90)
plt.xlabel("Tissue to train EMS", size=16)
plt.ylabel("Enrichment of EMS in Whole Blood", size=14)
plt.xlim([-0.5,49+0.5])
plt.ylim([0,4.5])
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/fig2_enr_wb.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()

#then, all tissue vs all tissue
gtex_colors.index = gtex_colors.tissue_abbrv
cnt = 0
vs = []
orders = []
aves= []
meds = []
colors = gtex_colors.loc[enr_med.index,"tissue_color_hex"]
for i in enr_med.index:
    l = pd.DataFrame(enr_med.loc[i,:].sort_values(ascending=False))
    ave = l.mean().values[0]
    med = l.median().values[0]
    v = l.loc[i,i]
    best = l.loc[i,:].sort_values(ascending=False).head(1).values[0]
    best_tissue = l.loc[i,:].sort_values(ascending=False).index[0]
    l["idx"] = np.arange(l.shape[0])+1
    order = l[l.index==i].idx.values[0]
    #color = gtex_colors.loc[l.index,"tissue_color_hex"]#この情報はshapeでencodeする
    vs.append(v)
    orders.append(order)
    aves.append(ave)
    meds.append(med)
 df = pd.DataFrame({"v":vs, "ave":aves, "color":colors, 
                   "name":list(enr_med.index), "order":orders, 
                  "med":meds}).sort_values(by="v", ascending=False)
#plot:
x = np.arange(df.shape[0])
y = df.v
plt.figure(figsize=(10, 6))
plt.scatter(x, y, color = df.color)
plt.scatter(x, df.ave, color = "black", marker="D", label="average across tissues")
tks = ("(rk=" + df["order"].astype(str) + ") " + df.index).str.replace('\\(rk\\=1\\)',"")
plt.xticks(x, tks, rotation=90)
plt.xlabel("Tissue to train EMS\n(rank, if not the best)", size=18)
plt.ylabel("Enrichment of EMS\nin matched tissue", size=18)
plt.xlim([-0.5,49+0.5])
plt.legend()
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/fig2_enr_all.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=400)
plt.show()
