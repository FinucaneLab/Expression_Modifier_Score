# -*- coding: utf-8 -*-
__author__ = 'QingboWang'
import sys
import pandas as pd
import numpy as np
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')

tissue_name = sys.argv[1]
leftout_chr = sys.argv[2] #with the name "chr" on


fn = "gs://qingbowang/ems_v1_test/{0}_negative_training_vg_annotated.tsv".format(tissue_name)
pos = hl.import_table(fn, impute=True)
pos = pos.filter(pos.vg.split("_")[0]!=leftout_chr)#filter to other chrs
pos = pos.to_pandas()
pos.index = pos.vg
del pos["vg"]
fn = "gs://qingbowang/ems_v1_test/{0}_negative_training_vg_annotated.tsv".format(tissue_name)
neg = hl.import_table(fn, impute=True)
neg = neg.filter(neg.vg.split("_")[0]!=leftout_chr)#filter to other chrs
neg = neg.to_pandas()
neg.index = neg.vg
del neg["vg"]




#remove the binary features for now
bin_features = ['DHS_Trynka', 'UTR_3_UCSC', 'UTR_5_UCSC', 'Coding_UCSC',
       'Intron_UCSC', 'TFBS_ENCODE', 'TSS_Hoffman', 'H3K27ac_PGC2',
       'H3K9ac_Trynka', 'Promoter_UCSC', 'H3K4me1_Trynka', 'H3K4me3_Trynka',
       'Enhancer_Hoffman', 'Transcribed_Hoffman']
pos.drop(bin_features, axis=1, inplace=True)
neg.drop(bin_features, axis=1, inplace=True)

def abs_log1p(x):
    return (np.log2(abs(x) + 1))
pos = pos.apply(lambda x: abs_log1p(x))
neg = neg.apply(lambda x: abs_log1p(x))

#and then ML -- no split since we are just looking at the feature importances
from sklearn.ensemble import RandomForestClassifier
y = [1]*pos.shape[0] + [0]*neg.shape[0]
X = pd.concat([pos, neg], axis=0)
regr = RandomForestClassifier(n_estimators=100, min_samples_split=2, max_features=0.5, random_state=1,
                              n_jobs=-1) #start with 100 to see how it goes, even tho definitely not enough
regr.fit(X, y)
importances = regr.feature_importances_
imp = pd.Series(importances)
imp.index = X.columns
imp.sort_values(inplace=True, ascending=False)

fn = "gs://qingbowang/ems_v1_test/bas_feat_imp_100trees_{0}_{1}_leftout.tsv".format(tissue_name, leftout_chr)
with hl.hadoop_open(fn, 'w') as f:
    imp.to_csv(f, sep="\t")