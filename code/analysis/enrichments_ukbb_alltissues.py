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

#aggregation:
fn = "gs://qingbowang/ems_v1_test/ukbb_trait_categ_mapper.tsv"
with hl.hadoop_open(fn, 'r') as f:
    mapper = pd.read_csv(f, sep="\t") #trait mapper
allcateg = list(mapper.categ.unique())
mu0 = {}
mu1 = {}
sd = {}
n = {}
for tissue_name in tissues[:2]: #unit test
    mu0[tissue_name] = []
    mu1[tissue_name] = []
    sd[tissue_name] = []
    n[tissue_name] = []
    ems = hl.read_table("gs://qingbowang/ems_v1_test/ems_pcausal_gtexvg_all{0}.ht".format(tissue_name))
    ems = ems.annotate(hg38_ID = ems.vg.split("_")[0] + "_" + ems.vg.split("_")[1] + "_" +ems.vg.split("_")[2] + "_" +ems.vg.split("_")[3]).key_by("hg38_ID").select("p_causal", "confidence_gain")
    for categ in allcateg:
        comp = hl.read_table("gs://qingbowang/UKBB_nc_pp_susie_maxpip_{0}.ht".format(categ))
        ht0 = ems.join(comp, how="left")
        ht0 = ht0.filter(hl.is_defined(ht0.max_pip)) #remove those that does not contain complex trait information
        mu0[tissue_name].append(ht0.aggregate(hl.agg.mean(ht0.max_pip)))
        ht = ht0.filter(ht0.confidence_gain>10)#filter by confidence gain here to reduce n for order_by
        if ht.count()<10000:
            ht = ht0.filter(ht0.confidence_gain>1)
        if ht.count()<10000:
            ht = ht0.filter(ht0.confidence_gain>0.1)
        ht = ht.order_by(hl.desc(ht.p_causal))
        ht = ht.add_index()
        ht = ht.annotate(top10k = ht.idx<10000)
        ht = ht.filter(ht.top10k)
        st = ht.aggregate(hl.agg.stats(ht.max_pip))
        mu1[tissue_name].append(st.mean)
        sd[tissue_name].append(st.stdev)
        n[tissue_name].append(st.n)
        print ("done {0}, {1}".format(categ, tissue_name))
        print (mu0[tissue_name])
        print (st)
mu0 = pd.DataFrame(mu0)
mu1 = pd.DataFrame(mu1)
sd = pd.DataFrame(sd)
n = pd.DataFrame(n)
mu0.index = allcateg
mu1.index = allcateg
sd.index = allcateg
n.index = allcateg
fn = "gs://qingbowang/comp_trait_roc/allpairs_mean_comp_pip_all.tsv"
with hl.hadoop_open(fn, 'w') as f:
    mu0.to_csv(f, sep="\t")
fn = "gs://qingbowang/comp_trait_roc/allpairs_mean_comp_pip_top10k.tsv"
with hl.hadoop_open(fn, 'w') as f:
    mu1.to_csv(f, sep="\t")
fn = "gs://qingbowang/comp_trait_roc/allpairs_sd_comp_pip_top10k.tsv"
with hl.hadoop_open(fn, 'w') as f:
    sd.to_csv(f, sep="\t")
fn = "gs://qingbowang/comp_trait_roc/allpairs_n_comp_pip_top10k.tsv"
with hl.hadoop_open(fn, 'w') as f:
    n.to_csv(f, sep="\t")


#plotting:
mu0 = []
mu1 = []
sd = []
n = []
for tissue_name in tissues:
    try:
        fn = "gs://qingbowang/comp_trait_roc/allpairs_mean_comp_pip_all_{0}.tsv".format(tissue_name)
        with hl.hadoop_open(fn, 'r') as f:
            mu0.append(pd.read_csv(f, sep="\t", index_col=0))
        fn = "gs://qingbowang/comp_trait_roc/allpairs_mean_comp_pip_top10k_{0}.tsv".format(tissue_name)
        with hl.hadoop_open(fn, 'r') as f:
            mu1.append(pd.read_csv(f, sep="\t", index_col=0))
        fn = "gs://qingbowang/comp_trait_roc/allpairs_sd_comp_pip_top10k_{0}.tsv".format(tissue_name)
        with hl.hadoop_open(fn, 'r') as f:
            sd.append(pd.read_csv(f, sep="\t", index_col=0))
        fn = "gs://qingbowang/comp_trait_roc/allpairs_n_comp_pip_top10k_{0}.tsv".format(tissue_name)
        with hl.hadoop_open(fn, 'r') as f:
            n.append(pd.read_csv(f, sep="\t", index_col=0))
    except:
        pass
mu0 = pd.concat(mu0, axis=1)
mu1 = pd.concat(mu1, axis=1)
sd = pd.concat(sd, axis=1)
n = pd.concat(n, axis=1)
enr = mu1 / mu0
err = sd / n.applymap(lambda x: np.sqrt(x)) / mu0


#color
fn = "gs://qingbowang/gtex_colors.tsv"
with hl.hadoop_open(fn, 'r') as f:
        gtex_colors = pd.read_csv(f, sep="\t", index_col=0)
gtex_colors["tissue_color_hex"] = "#" + gtex_colors["tissue_color_hex"].astype(str) 
fn = "gs://qingbowang/gtex_name_corres.csv"
with hl.hadoop_open(fn, 'r') as f:
        gtex_names = pd.read_csv(f, sep=",", index_col=0)
gtex_colors = pd.concat([gtex_colors, gtex_names], axis=1).fillna("0")
gtex_colors.index = gtex_colors.my_name.str.replace(" ", "")

colors = gtex_colors.loc[enr.columns,"tissue_color_hex"]
abb = gtex_colors.loc[enr.columns,"tissue_abbrv"]
#n_pos in random setting
fn = "gs://qingbowang/comp_trait_roc/n_pos_and_neg.tsv"
with hl.hadoop_open(fn, 'r') as f:
    n_all = pd.read_csv(f, sep="\t", index_col=0)    

#plot the enrichment
fig = plt.figure(figsize=(15,12))
ax = {}
i = 1
for trait in n_all.loc[:, "pos"].sort_values(ascending=False).index:
    ax[i] = fig.add_subplot(4, 3, i)
    lr = enr.loc[trait,:].sort_values()
    l = lr.tail(5) #should have named l and r opposite, but fine..
    r = lr.head(5)
    lerr = err.loc[trait,l.index]
    rerr = err.loc[trait,r.index]
    xtk_l = abb[lr.index].tail(5)
    xtk_r = abb[lr.index].head(5)
    cl = colors[lr.index].tail(5)
    cr = colors[lr.index].head(5)
    ax[i].set_title(trait)
    ax[i].errorbar(np.arange(5), list(r), list(rerr), ecolor=cr, fmt="none")
    ax[i].scatter(np.arange(5), list(r), color=cr)
    ax[i].errorbar(np.arange(7,12), list(l), list(lerr), ecolor=cl, fmt="none")
    ax[i].scatter(np.arange(7,12), list(l), color=cl)
    ax[i].set_xticks(list(np.arange(5))+list(np.arange(7,12)))
    ax[i].set_xticklabels(list(xtk_r)+list(xtk_l), rotation=60)
    i += 1
fig.text(0.5, 0, 'Tissue to train EMS', ha='center', size=22)
fig.text(0, 0.5, 'Trait PIP enrichment', va='center', rotation='vertical', size=22)
    
plt.tight_layout(pad= 2.0)
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/fig4_pipratio_allvsall_selected.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=250)

plt.show()


