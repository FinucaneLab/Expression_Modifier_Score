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

#plot x=random forest output, y=ems
tissue_name = "Whole_Blood"
fn = "gs://qingbowang/ems_v1_test/ems_p_causal_{0}.tsv".format(tissue_name)
with hl.hadoop_open(fn, 'r') as f:
    df = pd.read_csv(f, sep="\t", index_col=0)
plt.rcParams.update({'font.size': 18})
plt.figure(figsize=(6, 4.5))
plt.plot(df.index, df.p_causal, color="black")
plt.xlabel("Random forest output score")
plt.ylabel("EMS")
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/ems_vs_rf_output_score.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=500)
plt.show()    

#plot F1 measures for ROADMAP features
with hl.hadoop_open("gs://qingbowang/ems_v1_test/{0}_roadmapannot_f1measures.tsv".format(tissue_name), 'r') as f:
    df = pd.read_csv(f, sep="\t", index_col=0)
    plt.figure(figsize=(6, 4.5))
df.sort_values(by="f1", inplace=True, ascending=False)
plt.scatter(np.arange(df.shape[0]), df.f1, color="black")
plt.axvline(x=df.shape[0]*0.025, linestyle="--", color="black")
plt.xlabel("Feature")
plt.ylabel("F1 measure")
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/rm_f1measure.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=500)
    
#plot feature importances for basenji features
fn = "gs://qingbowang/ems_v1_test/bas_feat_imp_100trees_Whole_Blood.tsv"
with hl.hadoop_open(fn, 'r') as f:
    bs = pd.read_csv(f, sep="\t", index_col=0) 
#(remove TSS distance from the features)
bs = bs.iloc[1:,:]
plt.figure(figsize=(6, 4.5))
plt.scatter(np.arange(bs.shape[0]), bs["0"], color="black")
plt.axvline(x=100, linestyle="--", color="black")
plt.xlabel("Feature")
plt.ylabel("MDI")
plt.tight_layout()
fn = "gs://qingbowang/ems_v1_test/agg_for_fig1/bas_featimp.png"
with hl.hadoop_open(fn, 'wb') as f:
    plt.savefig(f, dpi=500)

    
#plot EMS per tissue:
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
            "Vagina"]#all 49 tissues
rocs = []
prcs = []
for tissue_name in tissues:
    fn = "gs://qingbowang/ems_v1_test/grid_search_result_{0}.tsv".format(tissue_name)
    with hl.hadoop_open(fn, 'r') as f: 
        df = pd.read_csv(f, sep="\t")
    rocs.append(df.auroc.values[0])
    prcs.append(df.auprc.values[0])
roc = pd.DataFrame({"roc":rocs,"prc":prcs})
roc.index = tissues
#add the sample size...
npos = []
for tissue_name in tissues:
    fn = "gs://qingbowang/ems_v1_test/{0}_positive_training_vg_beforeannot.tsv".format(tissue_name)
    with hl.hadoop_open(fn, 'r') as f:
            df = pd.read_csv(f, sep="\t")
    npos.append(df.shape[0])
roc["n_pos"] = npos
#plot:
from scipy import stats
x = roc.n_pos
y = roc.roc
plt.scatter(x,y)
(r,p) = stats.pearsonr(x, y)
plt.xlabel("n(positive sample)")
plt.ylabel("AUROC")
plt.title("(pearson r, p-value) = ({0},{1})".format(str(r)[:6], str(p)[:6]))
plt.show()
