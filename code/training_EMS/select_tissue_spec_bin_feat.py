# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

import sys
import numpy as np
import pandas as pd
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')

tissue_name = sys.argv[1]
marks = ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3" ,"H3K9ac", "H3K9me3"]

df_prec = []
df_rec = []
df = []
for m in marks:
    pos = hl.import_table("gs://qingbowang/ems_v1_test/{0}_positive_training_vg_roadmapannot_{1}.tsv".format(tissue_name, m), impute=True).to_pandas()
    neg = hl.import_table("gs://qingbowang/ems_v1_test/{0}_negative_training_vg_roadmapannot_{1}.tsv".format(tissue_name, m), impute=True).to_pandas()
    #set the index_col
    pos.index = pos.vg
    neg.index = neg.vg
    del pos["vg"]
    del neg["vg"]
    #for each of the column, get the true vs false, and that easily gives us f1 measure -> then compare
    TP = pos.sum(axis=0)
    FN = pos.shape[0] - TP
    FP = neg.sum(axis=0)
    TN = neg.shape[0] - FP #not used
    #and downsample the negative to fit the 1:1 ratio
    FP_weight = FP * (pos.shape[0]/neg.shape[0])
    precision = TP / (TP + FP_weight)
    recall = TP / (TP + FN)
    dfsub = pd.concat([TP, FN, FP, TN, precision, recall], axis=1)
    df.append(dfsub)
df = pd.concat(df, axis=0)
df.columns = ["TP", "FN", "FP", "TN", "precision", "recall"]
df["f1"] = 2/(df.precision**-1 + df.recall**-1)
df["OR"] = (df.TP / df.FN) / (df.FP / df.TN) #= how much it is enriched in the positive set

#and save
with hl.hadoop_open("gs://qingbowang/ems_v1_test/{0}_roadmapannot_f1measures.tsv".format(tissue_name), 'w') as f:
    df.to_csv(f, sep="\t")

#tissue union binaries
df = []
cols_to_use_bin = ['DHS_Trynka', 'UTR_3_UCSC', 'UTR_5_UCSC', 'Coding_UCSC', 'Intron_UCSC', 'TFBS_ENCODE', 'TSS_Hoffman', 'H3K27ac_PGC2', 'H3K9ac_Trynka', 'Promoter_UCSC', 'H3K4me1_Trynka', 'H3K4me3_Trynka', 'Enhancer_Hoffman', 'Transcribed_Hoffman']
pos = hl.import_table("gs://qingbowang/ems_v1_test/{0}_positive_training_vg_beforeannot.tsv".format(tissue_name), impute=True)
pos = pos.select(*cols_to_use_bin).to_pandas()
neg = hl.import_table("gs://qingbowang/ems_v1_test/{0}_negative_training_vg_beforeannot.tsv".format(tissue_name), impute=True)
neg = neg.select(*cols_to_use_bin).to_pandas()
TP = pos.sum(axis=0)
FN = pos.shape[0] - TP
FP = neg.sum(axis=0)
TN = neg.shape[0] - FP #not used
#and downsample the negative to fit the 1:1 ratio
FP_weight = FP * (pos.shape[0]/neg.shape[0])
precision = TP / (TP + FP_weight)
recall = TP / (TP + FN)
dfsub = pd.concat([TP, FN, FP, TN, precision, recall], axis=1)
df.append(dfsub)
df = pd.concat(df, axis=0)
df.columns = ["TP", "FN", "FP", "TN", "precision", "recall"]
df["f1"] = 2/(df.precision**-1 + df.recall**-1)
df["OR"] = (df.TP / df.FN) / (df.FP / df.TN) #= how much it is enriched in the positive set

#and save
with hl.hadoop_open("gs://qingbowang/ems_v1_test/{0}_tissue_union_annot_f1measures.tsv".format(tissue_name), 'w') as f:
    df.to_csv(f, sep="\t")
