# -*- coding: utf-8 -*-
__author__ = 'QingboWang'
import numpy as np
import pandas as pd
import sys

tissue_name = sys.argv[1]
pos = pd.read_csv("{0}_positive_training_vg_annotated.tsv".format(tissue_name), sep="\t", index_col=0)
neg = pd.read_csv("{0}_negative_training_vg_annotated.tsv".format(tissue_name), sep="\t", index_col=0)

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

#and then feature selection in RF -- no split since we are just looking at the feature importances
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
imp.to_csv("bas_feature_imp_100forest_{0}.tsv".format(tissue_name), sep="\t")