# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

import sys
import pandas as pd
import numpy as np
import time as tm
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics

def abs_log1p(x):
    return (np.log2(abs(x) + 1))

tissue_name = sys.argv[1]
leftout_chr = sys.argv[2]

#subset to features to use
basfeats = pd.read_csv("bas_feature_top100_touse_per_leftout_chr.tsv", sep="\t", index_col=0)
rmfeats = pd.read_csv("roadmapannot_feat_to_use_per_leftout_chr.tsv", sep="\t", index_col=0)
feat_to_use_bas = basfeats.loc[:,leftout_chr]
feat_to_use_rm = rmfeats.loc[:,leftout_chr].dropna()
feat_to_use_bin = ['DHS_Trynka', 'UTR_3_UCSC', 'UTR_5_UCSC', 'Coding_UCSC',
       'Intron_UCSC', 'TFBS_ENCODE', 'Promoter_UCSC', 'Enhancer_Hoffman'] #all the non-histone features

#load data
pos = []
neg = []
posbas = pd.read_csv("{0}_positive_training_vg_annotated.tsv".format(tissue_name), sep="\t", index_col=0)
negbas = pd.read_csv("{0}_negative_training_vg_annotated.tsv".format(tissue_name), sep="\t", index_col=0)

#filter to leave one chr
posbas = posbas[posbas.index.str.split("_").str[0]!=leftout_chr]
negbas = negbas[negbas.index.str.split("_").str[0]!=leftout_chr]

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

    #filter to leave one chr
    posrm = posrm[posrm.index.str.split("_").str[0]!=leftout_chr]
    negrm = negrm[negrm.index.str.split("_").str[0]!=leftout_chr]

    posrm = posrm.loc[:, feat_to_use_rm[feat_to_use_rm.str.contains(m)]]
    negrm = negrm.loc[:, feat_to_use_rm[feat_to_use_rm.str.contains(m)]]
    posrm = posrm*1
    negrm = negrm*1
    pos.append(posrm)
    neg.append(negrm)
pos = pd.concat(pos, axis=1)
neg = pd.concat(neg, axis=1)

#re-order, just for convenience -- no need
#feat_to_use = feat_to_use_bin + list(feat_to_use_rm) + list(feat_to_use_bas) #結論から言うとこれをちゃんとやっておけばよかった....
#pos = pos.loc[:,feat_to_use]

feat_to_use = list(pos.columns)

#define the performance test

def cv_rf_aucs(pos_X, neg_X, feat_to_use, nest=100, min_split=2, max_feat=0.1, max_depth=None, criterion="gini", fold=10, frac_train=0.3, seed=1):
    #returns the AUPRC and AUROC under 1:1 ratio, average over the folds
    print ("starting param: nest={0}, minsplit={1}, maxfeat={2}, {3}".format(nest, min_split, max_feat, tm.ctime()))
    np.random.seed(seed)
    pos_X = pos_X.loc[:, feat_to_use]
    neg_X = neg_X.loc[:, feat_to_use]
    aurocs = []
    auprcs = []
    for iter in range(fold):
        train_pos = np.random.choice(2, pos_X.shape[0], p=[frac_train, 1-frac_train])
        train_neg = np.random.choice(2, neg_X.shape[0], p=[frac_train, 1-frac_train])
        X_train = pd.concat([pos_X[train_pos == 1], neg_X[train_neg == 1]])
        X_test = pd.concat([pos_X[train_pos == 0], neg_X[train_neg == 0]])
        y_train = np.array([1] * sum(train_pos == 1) + [0] * sum(train_neg == 1))
        y_test = np.array([1] * sum(train_pos == 0) + [0] * sum(train_neg == 0))
        regr = RandomForestClassifier(n_estimators=nest, min_samples_split=min_split,  max_features=max_feat, random_state=1, criterion=criterion, n_jobs=-1)
        if len(feat_to_use)==1: #need to reshape when only one features
            regr.fit(np.array(X_train).reshape(-1, 1), y_train,
                     sample_weight=np.array([1] * sum(train_pos == 1) + [sum(train_pos == 1) / sum(train_neg == 1)] * sum(train_neg == 1)))
        else:
            regr.fit(X_train, y_train,
                     sample_weight=np.array([1] * sum(train_pos == 1) + [sum(train_pos == 1) / sum(train_neg == 1)] * sum(train_neg == 1)))
        if len(feat_to_use) == 1: #need to reshape when only one features
            y_prob = regr.predict_proba(np.array(X_test).reshape(-1, 1))
        else:
            y_prob = regr.predict_proba(X_test)
        y_prob = y_prob[:, 1]
        #get the auroc
        fpr, tpr, threshold = metrics.roc_curve(y_test, y_prob,
                               sample_weight=np.array([1] * sum(train_pos == 0) + [sum(train_pos==0)/sum(train_neg==0)] * sum(train_neg == 0)))
        auroc = metrics.auc(fpr, tpr)
        aurocs.append(auroc)
        #and auprc
        auprc = metrics.average_precision_score(y_test, y_prob,
                                                sample_weight=np.array([1] * sum(train_pos == 0) + [sum(train_pos==0)/sum(train_neg==0)] * sum(train_neg == 0)))
        auprcs.append(auprc)
        print ("done {0}th iter, {1}".format(iter, tm.ctime()))
        if iter<2:
            print ("auprc: {0}".format(auprc)) #sanity check
    #get the best and record that, as well as update the roc
    aurocs_ave = np.array(aurocs).mean()
    auprcs_ave = np.array(auprcs).mean()
    return ([aurocs_ave, auprcs_ave])


def random_search(seed=2020, size=100):
    np.random.seed(seed)
    rocs = []
    prcs = []
    #nests = np.random.choice([10, 32, 100, 316, 1000, 3162], size=size)
    nests = np.random.choice([10, 32, 100, 316, 1000], size=size) #max 1000, since it takes too long
    msps = np.random.choice([2, 8, 16], size=size)
    mfeats = np.random.choice([0.1, 0.32, 0.64, 1], size=size)
    mdepths = np.random.choice([None, 2], size=size)
    criterias = np.random.choice(["gini", "entropy"], size)
    for i in range(size):
        roc, prc = cv_rf_aucs(pos, neg, feat_to_use, nest=nests[i], min_split=msps[i], max_feat=mfeats[i],
                                              max_depth=mdepths[i], criterion=criterias[i])
        rocs.append(roc)
        prcs.append(prc)
        print ("done {0} th random search out of {1}, {2}".format(i, size, tm.ctime()))
    df = pd.DataFrame({"auroc": rocs, "auprc": prcs, "nest": nests,
                       "msp": msps, "mfeat": mfeats, "mdepth": mdepths, "crit":criterias})
    return (df)

df = random_search(size=100)
auroc, auprc, nest, msp, mfeat, mdepth, criteria = list(df.sort_values(by="auroc", ascending=False).iloc[0,:])

print ("best:")
print (df.sort_values(by="auroc", ascending=False).head(1))

#then the grid search:

nests_grid = (np.arange(0.5, 1.6, .1) * nest).astype(int)
msp_grid = np.arange(max(2,msp/2), msp*1.5+0.1, 2).astype(int)
mfeat_grid = np.arange(0.5, min(2/mfeat, 1)+0.01, 0.1) * mfeat

def grid_search(nests_grid, msps_grid, mfeats_grid, criterion="Gini", mdepth="None"):
    rocs = []
    prcs = []
    nests = []
    msps = []
    mfeats = []
    for nest in nests_grid:
        for msp in msps_grid:
            for mfeat in mfeats_grid:
                roc, prc = cv_rf_aucs(pos, neg, feat_to_use, nest=nest, min_split=msp, max_feat=mfeat,
                                              max_depth=mdepth, criterion=criterion)
                rocs.append(roc)
                prcs.append(prc)
                nests.append(nest)
                msps.append(msp)
                mfeats.append(mfeat)
                print ("done nest={0}, msp={1}, mfeat={2}, {3}".format(nest, msp, mfeat, tm.ctime()))
    df = pd.DataFrame({"auroc": rocs, "auprc": prcs, "nest": nests, "msp": msps,
                       "mfeat": mfeats, "mdepth_fixed":mdepth, "criteria_fixed":criteria})
    return (df)

print ("starting grid search")
dfgrid = grid_search(nests_grid, msp_grid, mfeat_grid, criterion=criteria, mdepth=mdepth)


dfgrid.sort_values(by="auroc", ascending=False).to_csv("grid_search_result_{0}_{1}_left_out.tsv".format(tissue_name, leftout_chr), sep="\t")
df.sort_values(by="auroc", ascending=False).to_csv("random_search_result_{0}_{1}_left_out.tsv".format(tissue_name, leftout_chr), sep="\t")
