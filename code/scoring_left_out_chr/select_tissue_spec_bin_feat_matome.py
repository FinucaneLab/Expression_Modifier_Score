# -*- coding: utf-8 -*-
__author__ = 'QingboWang'


import sys
import numpy as np
import pandas as pd
from scipy import stats
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')


chr_list = []
for i in range(1,22+1):
    chr_list.append("chr" + str(i))
chr_list.append("chrX")
tissue_name = "Whole_Blood"
#and draw some figures
from matplotlib import pyplot as plt
from matplotlib import cm
def fisher_p(l):
    OR, pval = stats.fisher_exact([[l["TP"], l["FN"]], [l["FP"], l["TN"]]])
    return (pval)
marks = ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3" ,"H3K9ac", "H3K9me3"]
features_to_use_final = {}
for chr in chr_list: #leftout chr
    with hl.hadoop_open("gs://qingbowang/ems_v1_test/{0}_roadmapannot_f1measures_{1}_leftout.tsv".format(tissue_name, chr), 'r') as f:
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
    features_to_use_final[chr] = features_to_use

    #export the meta data
    with hl.hadoop_open("gs://qingbowang/ems_v1_test/{0}_roadmapannot_selectedfeatures_summary_{1}.tsv".format(tissue_name, chr), 'w') as f:
        df.loc[features_to_use,:].to_csv(f, sep="\t")

features_to_use_final = pd.DataFrame.from_dict(features_to_use_final, orient='index').fillna("NA").T
with hl.hadoop_open("gs://qingbowang/ems_v1_test/roadmapannot_feat_to_use_per_leftout_chr.tsv", 'w') as f:
    features_to_use_final.to_csv(f, sep="\t")




