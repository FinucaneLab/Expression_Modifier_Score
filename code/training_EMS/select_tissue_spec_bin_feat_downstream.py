# -*- coding: utf-8 -*-
__author__ = 'QingboWang'


import sys
import numpy as np
import pandas as pd
from scipy import stats
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')
tissue_list = ["Whole_Blood",
                "Muscle_Skeletal",
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
                "Liver",
                "Adipose_Visceral_Omentum",
                "Pancreas",
                "Lung",
                "Pituitary",
                "Brain_Nucleus_accumbens_basal_ganglia",
                "Colon_Transverse",
                "Adipose_Subcutaneous",
                "Esophagus_Mucosa",
                "Brain_Cerebellum",
                "Brain_Cortex",
                "Thyroid",
                "Stomach",
                "Breast_Mammary_Tissue",
                "Colon_Sigmoid",
                "Skin_Not_Sun_Exposed_Suprapubic",
                "Testis"]
#and draw some figures
from matplotlib import pyplot as plt
from matplotlib import cm
def fisher_p(l):
    OR, pval = stats.fisher_exact([[l["TP"], l["FN"]], [l["FP"], l["TN"]]])
    return (pval)
marks = ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3" ,"H3K9ac", "H3K9me3"]
features_to_use_final = {}
for tissue_name in tissue_list:
    with hl.hadoop_open("gs://qingbowang/ems_v1_test/{0}_roadmapannot_f1measures.tsv".format(tissue_name), 'r') as f:
        df = pd.read_csv(f, sep="\t", index_col=0)
    df.sort_values(by="f1", ascending=False, inplace=True)
    # choose the top 2.5% of f1
    top_f1 = df.iloc[:int(np.ceil(df.shape[0]*0.025)),:].index
    # choose the top 2.5% of precision
    top_prec = df.sort_values(by="precision", ascending=False).iloc[:int(np.ceil(df.shape[0]*0.025)),:].index
    # the union of those will be the ones to select
    tops = np.union1d(top_f1, top_prec)
    # choose bottom 5% of f1
    df_bottom_f1 = df.iloc[-int(np.floor(df.shape[0]*0.025)):,]
    # filter by OR<0.5, pval<0.05
    df_bottom_f1 = df_bottom_f1[df_bottom_f1.OR<0.5]
    df_bottom_f1["fisher_p"] = df_bottom_f1.apply(lambda x: fisher_p(x), axis=1)
    df_bottom_f1 = df_bottom_f1[df_bottom_f1.fisher_p<0.05]
    # those will be the ones to select
    bottoms = df_bottom_f1.index

    features_to_use = np.union1d(tops, bottoms)
    #save this for final export
    features_to_use_final[tissue_name] = features_to_use

    #export the meta data
    with hl.hadoop_open("gs://qingbowang/ems_v1_test/{0}_roadmapannot_selectedfeatures_summary.tsv".format(tissue_name), 'w') as f:
        df.loc[features_to_use,:].to_csv(f, sep="\t")
    #and plot the results for later tracking
    i = 0
    for m in marks:
        dfsub = df[df.index.str.contains(m)]
        dfsub["idx"] = np.arange(dfsub.shape[0])+i
        plt.scatter(dfsub.idx, dfsub.f1, color=cm.tab10(i), label=m)
        i += 1
    plt.legend()
    plt.title(tissue_name)
    plt.xlabel("features, sorted by f1 score")
    plt.ylabel("f1 score")
    with hl.hadoop_open("gs://qingbowang/ems_v1_test/figures/{0}_f1_per_histonmark.png".format(tissue_name), 'wb') as f:
        plt.savefig(f, dpi=400)
    #precision
    i = 0
    dfp = df.sort_values(by="precision", ascending=False)
    for m in marks:
        dfsub = dfp[dfp.index.str.contains(m)]
        dfsub["idx"] = np.arange(dfsub.shape[0])+i
        plt.scatter(dfsub.idx, dfsub.precision, color=cm.tab10(i), label=m)
        i += 1
    plt.legend()
    plt.title(tissue_name)
    plt.xlabel("features, sorted by precision")
    plt.ylabel("precision")
    with hl.hadoop_open("gs://qingbowang/ems_v1_test/figures/{0}_precision_per_histonmark.png".format(tissue_name), 'wb') as f:
        plt.savefig(f, dpi=400)

    #overall
    plt.scatter(np.arange(df.shape[0])/df.shape[0], df.f1, color="black", s=8)
    plt.vlines(x=0.025, ymin=0, ymax=1)
    plt.title(tissue_name)
    plt.xlabel("features, sorted by f1 score (n={0})".format(df.shape[0]))
    plt.ylabel("f1 score")
    with hl.hadoop_open("gs://qingbowang/ems_v1_test/figures/{0}_f1_overall.png".format(tissue_name), 'wb') as f:
        plt.savefig(f, dpi=400)

    plt.scatter(np.arange(dfp.shape[0])/dfp.shape[0], dfp.precision, color="black", s=8)
    plt.vlines(x=0.025, ymin=0, ymax=1)
    plt.xlabel("features, sorted by precision (n={0})".format(df.shape[0]))
    plt.ylabel("precision")
    with hl.hadoop_open("gs://qingbowang/ems_v1_test/figures/{0}_precision_overall.png".format(tissue_name), 'wb') as f:
        plt.savefig(f, dpi=400)
    
    print ("done {0}".format(tissue_name))


features_to_use_final = pd.DataFrame.from_dict(features_to_use_final, orient='index').fillna("NA").T
with hl.hadoop_open("gs://qingbowang/ems_v1_test/roadmapannot_feat_to_use_all.tsv", 'w') as f:
    features_to_use_final.to_csv(f, sep="\t")







