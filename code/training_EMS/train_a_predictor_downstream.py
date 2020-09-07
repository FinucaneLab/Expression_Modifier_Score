# -*- coding: utf-8 -*-
__author__ = 'QingboWang'


import sys
import pandas as pd
import numpy as np
import time as tm
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import subprocess
from matplotlib import pyplot as plt

#load the data (load from the cloud first)
tissue_name = sys.argv[1]

subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/grid_search_result_{0}.tsv".format(tissue_name), "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_positive_training_vg_annotated.tsv".format(tissue_name), "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_negative_training_vg_annotated.tsv".format(tissue_name), "./"])

subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_negative_training_vg_roadmapannot_H3K27ac.tsv".format(tissue_name), "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_negative_training_vg_roadmapannot_H3K27me3.tsv".format(tissue_name), "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_negative_training_vg_roadmapannot_H3K36me3.tsv".format(tissue_name), "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_negative_training_vg_roadmapannot_H3K4me1.tsv".format(tissue_name), "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_negative_training_vg_roadmapannot_H3K4me3.tsv".format(tissue_name), "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_negative_training_vg_roadmapannot_H3K9ac.tsv".format(tissue_name), "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_negative_training_vg_roadmapannot_H3K9me3.tsv".format(tissue_name), "./"])

subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_positive_training_vg_roadmapannot_H3K27ac.tsv".format(tissue_name), "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_positive_training_vg_roadmapannot_H3K27me3.tsv".format(tissue_name), "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_positive_training_vg_roadmapannot_H3K36me3.tsv".format(tissue_name), "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_positive_training_vg_roadmapannot_H3K4me1.tsv".format(tissue_name), "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_positive_training_vg_roadmapannot_H3K4me3.tsv".format(tissue_name), "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_positive_training_vg_roadmapannot_H3K9ac.tsv".format(tissue_name), "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/{0}_positive_training_vg_roadmapannot_H3K9me3.tsv".format(tissue_name), "./"])

subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/n_pos_and_neg_per_tissue.tsv", "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/bas_feature_top100_touse.tsv", "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/roadmapannot_feat_to_use_all.tsv", "./"])


def abs_log1p(x):
    return (np.log2(abs(x) + 1))


#subset to features to use
basfeats = pd.read_csv("bas_feature_top100_touse.tsv", sep="\t", index_col=0)
rmfeats = pd.read_csv("roadmapannot_feat_to_use_all.tsv", sep="\t", index_col=0)
feat_to_use_bas = basfeats.loc[:,tissue_name]
feat_to_use_rm = rmfeats.loc[:,tissue_name].dropna()
feat_to_use_bin = ['DHS_Trynka', 'UTR_3_UCSC', 'UTR_5_UCSC', 'Coding_UCSC',
       'Intron_UCSC', 'TFBS_ENCODE', 'Promoter_UCSC', 'Enhancer_Hoffman'] #all the non-histone features

#load data
pos = []
neg = []
posbas = pd.read_csv("{0}_positive_training_vg_annotated.tsv".format(tissue_name), sep="\t", index_col=0)
negbas = pd.read_csv("{0}_negative_training_vg_annotated.tsv".format(tissue_name), sep="\t", index_col=0)
posbas = posbas.loc[:,feat_to_use_bin + list(feat_to_use_bas)]
negbas = negbas.loc[:,feat_to_use_bin + list(feat_to_use_bas)]
#and turn the log1p
posbas = posbas.apply(lambda x: abs_log1p(x))
negbas = negbas.apply(lambda x: abs_log1p(x))

pos.append(posbas)
neg.append(negbas)
marks = ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3" ,"H3K9ac", "H3K9me3"]
for m in marks:
    posrm = pd.read_csv("{0}_positive_training_vg_roadmapannot_{1}.tsv".format(tissue_name, m), sep="\t", index_col=0)
    negrm = pd.read_csv("{0}_negative_training_vg_roadmapannot_{1}.tsv".format(tissue_name, m), sep="\t", index_col=0)
    posrm = posrm.loc[:, feat_to_use_rm[feat_to_use_rm.str.contains(m)]]
    negrm = negrm.loc[:, feat_to_use_rm[feat_to_use_rm.str.contains(m)]]
    posrm = posrm*1
    negrm = negrm*1
    pos.append(posrm)
    neg.append(negrm)
pos = pd.concat(pos, axis=1)
neg = pd.concat(neg, axis=1)

feat_to_use = list(pos.columns)


dfgrid = pd.read_csv("grid_search_result_{0}.tsv".format(tissue_name), sep="\t", index_col=0)
#remove the mdepth column, that was just by mistake
if "mdepth" in dfgrid.columns:
    del dfgrid["mdepth"]

dfgrid.fillna(1000000, inplace=True) #to fill NaN with basically None
print (dfgrid.head(2)) #sanity check


# train using 100% of the data, for the final EMS
auroc, auprc, nest, msp, mfeat, mdepth, criteria = list(dfgrid.sort_values(by="auroc", ascending=False).iloc[0,:]) #best parameter set.

regr = RandomForestClassifier(n_estimators=nest, min_samples_split=msp, max_features=mfeat, random_state=1,
                              max_depth=mdepth, criterion=criteria, n_jobs=-1)
X = pd.concat([pos, neg])
y = np.array([1] * pos.shape[0] + [0] * neg.shape[0])
regr.fit(X, y, sample_weight=np.array([1] * pos.shape[0] + [pos.shape[0] / neg.shape[0]] * neg.shape[0]))

#save
import pickle
filename = 'ems_rfmodel_{0}.sav'.format(tissue_name)
pickle.dump(regr, open(filename, 'wb'))


# do the scoring 100 times, 90% training + 10% test data

frac_test = 0.1 #10% for the test
fold = 100 #100 times to get a smooth dist.
y_real_all = []
y_prob_all = []
np.random.seed(1)
for i in range(fold):
    train_pos = np.random.choice(2, pos.shape[0], p=[frac_test, 1-frac_test])
    train_neg = np.random.choice(2, neg.shape[0], p=[frac_test, 1-frac_test])
    X_train = pd.concat([pos[train_pos == 1], neg[train_neg == 1]])
    X_test = pd.concat([pos[train_pos == 0], neg[train_neg == 0]])
    y_train = np.array([1] * sum(train_pos == 1) + [0] * sum(train_neg == 1))
    y_test = np.array([1] * sum(train_pos == 0) + [0] * sum(train_neg == 0))
    regr = RandomForestClassifier(n_estimators=nest, min_samples_split=msp, max_features=mfeat, random_state=1,
                              max_depth=mdepth, criterion=criteria, n_jobs=-1)
    #train the predictor
    regr.fit(X_train, y_train,
             sample_weight=np.array([1] * sum(train_pos == 1) + [sum(train_pos == 1) / sum(train_neg == 1)] * sum(train_neg == 1)))
    y_prob = regr.predict_proba(X_test)
    y_prob = y_prob[:,1]
    #save the prob
    y_real_all = y_real_all + list(y_test)
    y_prob_all = y_prob_all + list(y_prob)
    print ("{0} done {1}th iter out of {2}".format(tm.ctime(), i, fold))
#check distribution, and get the conversion (assuming that 90% ~= 100%)
pd.Series(np.array(y_prob_all)[np.array(y_real_all)==1]).hist(alpha=0.3, bins=100+1, color="r", density=True, label="positive")
pd.Series(np.array(y_prob_all)[np.array(y_real_all)==0]).hist(alpha=0.3, bins=100+1, color="b", density=True, label="negative")
plt.xlim([0,1])
plt.title("smoothed distribution after 100 iter, 90% training")
plt.legend()
plt.savefig("smooth_dist_{0}.png".format(tissue_name), dpi=400)
plt.clf()

#load the n_pos and n_others
ns = pd.read_csv("n_pos_and_neg_per_tissue.tsv", sep="\t", index_col=0)
n_pos = ns.loc[tissue_name,"n_pos"]
n_neg = ns.loc[tissue_name,"n_neg"] #actually having a threshold here seems fine.


f_pos = pd.Series(np.array(y_prob_all)[np.array(y_real_all)==1])\
            .value_counts().sort_index()/sum(np.array(y_real_all)==1)
f_neg = pd.Series(np.array(y_prob_all)[np.array(y_real_all)==0])\
            .value_counts().sort_index()/sum(np.array(y_real_all)==0)

#turn it to smoothed distribution
win = np.arange(0,1,0.001)
winpos = pd.Series(win).apply(lambda x: f_pos[(x-0.005<f_pos.index) & (f_pos.index<x+0.005)].sum())
winneg = pd.Series(win).apply(lambda x: f_neg[(x-0.005<f_neg.index) & (f_neg.index<x+0.005)].sum())

p_causal = (winpos*n_pos) / ( (winpos*n_pos) + (winneg*n_neg) )
p_causal = pd.DataFrame(p_causal)
p_causal.index = win
p_causal.columns = ["p_causal"]
p_causal["confidence_gain"] = p_causal.p_causal / (n_pos/n_neg)
p_causal["confidence_gain_log10"] = np.log10(p_causal.confidence_gain)

p_causal = p_causal[(p_causal.p_causal>0) & (p_causal.p_causal<1)] #pos neg 両方出てるやつじゃないとダメ. 信用できないので.

#save
p_causal.to_csv("ems_p_causal_{0}.tsv".format(tissue_name), sep="\t")

#copy to the bucket
subprocess.call(["gsutil", "cp" ,"./ems_rfmodel_{0}.sav".format(tissue_name) ,"gs://qingbowang/ems_v1_test/ems_rfmodel_{0}.sav".format(tissue_name)])
subprocess.call(["gsutil", "cp", "./ems_p_causal_{0}.tsv".format(tissue_name), "gs://qingbowang/ems_v1_test/ems_p_causal_{0}.tsv".format(tissue_name)])


#and remove the files:
subprocess.call(["rm", "./grid_search_result_{0}.tsv".format(tissue_name)])
subprocess.call(["rm", "./{0}_positive_training_vg_annotated.tsv".format(tissue_name)])
subprocess.call(["rm", "./{0}_negative_training_vg_annotated.tsv".format(tissue_name)])

subprocess.call(["rm", "./{0}_negative_training_vg_roadmapannot_H3K27ac.tsv".format(tissue_name)])
subprocess.call(["rm", "./{0}_negative_training_vg_roadmapannot_H3K27me3.tsv".format(tissue_name)])
subprocess.call(["rm", "./{0}_negative_training_vg_roadmapannot_H3K36me3.tsv".format(tissue_name)])
subprocess.call(["rm", "./{0}_negative_training_vg_roadmapannot_H3K4me1.tsv".format(tissue_name)])
subprocess.call(["rm", "./{0}_negative_training_vg_roadmapannot_H3K4me3.tsv".format(tissue_name)])
subprocess.call(["rm", "./{0}_negative_training_vg_roadmapannot_H3K9ac.tsv".format(tissue_name)])
subprocess.call(["rm", "./{0}_negative_training_vg_roadmapannot_H3K9me3.tsv".format(tissue_name)])

subprocess.call(["rm", "./{0}_positive_training_vg_roadmapannot_H3K27ac.tsv".format(tissue_name)])
subprocess.call(["rm", "./{0}_positive_training_vg_roadmapannot_H3K27me3.tsv".format(tissue_name)])
subprocess.call(["rm", "./{0}_positive_training_vg_roadmapannot_H3K36me3.tsv".format(tissue_name)])
subprocess.call(["rm", "./{0}_positive_training_vg_roadmapannot_H3K4me1.tsv".format(tissue_name)])
subprocess.call(["rm", "./{0}_positive_training_vg_roadmapannot_H3K4me3.tsv".format(tissue_name)])
subprocess.call(["rm", "./{0}_positive_training_vg_roadmapannot_H3K9ac.tsv".format(tissue_name)])
subprocess.call(["rm", "./{0}_positive_training_vg_roadmapannot_H3K9me3.tsv".format(tissue_name)])

