# -*- coding: utf-8 -*-
__author__ = 'QingboWang'
import numpy as np
import copy as cp
import sys, math, cmath
import time as tm
import pandas as pd
import os
import datetime
import random as rd
from matplotlib import pyplot as plt
from scipy import stats
import re
import mmap
import glob
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq


#official figure for the performance comparison of EMS
#to useだけに絞ってみるか....

#do this in the VM. it's going to be large..
"""in local
gcloud beta compute instances create "rf-vm" --machine-type=n1-highcpu-96 --scopes=storage-rw
gcloud beta compute --project "encode-uk-biobank" ssh --zone "us-central1-b" "rf-vm"
gsutil cp /Users/qingbowang/PycharmProjects/python3projects/gtex_finemapping/ems_pipe_test_20200125/20200820_rf_performance_comparison.py gs://qingbowang/20200820_rf_performance_comparison.py ./


"""
"""in vm
#and run this inside the ssh:
#load the package
sudo apt-get install bzip2 libxml2-dev
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
pip install --user sklearn
pip install --user pandas
#pip install --user pickle
pip install --user matplotlib

#copy the script from cloud
gsutil cp gs://qingbowang/20200820_rf_performance_comparison.py ./

#and other files to use
gsutil cp gs://qingbowang/ems_v1_test/bas_feature_top100_touse.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/roadmapannot_feat_to_use_all.tsv ./

tissue_name=Whole_Blood
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_annotated.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_annotated.tsv ./


gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_roadmapannot_H3K27ac.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_roadmapannot_H3K27me3.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_roadmapannot_H3K36me3.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_roadmapannot_H3K4me1.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_roadmapannot_H3K4me3.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_roadmapannot_H3K9ac.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_roadmapannot_H3K9me3.tsv ./

gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_roadmapannot_H3K27ac.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_roadmapannot_H3K27me3.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_roadmapannot_H3K36me3.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_roadmapannot_H3K4me1.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_roadmapannot_H3K4me3.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_roadmapannot_H3K9ac.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_roadmapannot_H3K9me3.tsv ./



"""


def abs_log1p(x):
    return (np.log2(abs(x) + 1))

tissue_name = "Whole_Blood"


#subset to features to use
basfeats = pd.read_csv("bas_feature_top100_touse.tsv", sep="\t", index_col=0)
rmfeats = pd.read_csv("roadmapannot_feat_to_use_all.tsv", sep="\t", index_col=0)
feat_to_use_bas = basfeats.loc[:,tissue_name][1:] #excludes tss
feat_tss = ["tss_distance"]
feat_to_use_rm = rmfeats.loc[:,tissue_name].dropna()
feat_to_use_bin = ['DHS_Trynka', 'UTR_3_UCSC', 'UTR_5_UCSC', 'Coding_UCSC',
       'Intron_UCSC', 'TFBS_ENCODE', 'Promoter_UCSC', 'Enhancer_Hoffman'] #all the non-histone features

#load data
pos = []
neg = []
posbas = pd.read_csv("{0}_positive_training_vg_annotated.tsv".format(tissue_name), sep="\t", index_col=0)
negbas = pd.read_csv("{0}_negative_training_vg_annotated.tsv".format(tissue_name), sep="\t", index_col=0)
posbas = posbas.loc[:,feat_to_use_bin + list(feat_to_use_bas) + feat_tss]
negbas = negbas.loc[:,feat_to_use_bin + list(feat_to_use_bas) + feat_tss] #let's use the full features
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
    negrm = negrm.loc[:, feat_to_use_rm[feat_to_use_rm.str.contains(m)]] #Let's use the full features
    posrm = posrm*1
    negrm = negrm*1
    pos.append(posrm)
    neg.append(negrm)
pos = pd.concat(pos, axis=1)
neg = pd.concat(neg, axis=1)

#fill na with 0 for simplicity
pos.fillna(0, inplace=True)
neg.fillna(0, inplace=True)
pos.replace([np.inf, -np.inf], 0, inplace=True)
neg.replace([np.inf, -np.inf], 0, inplace=True)
#and to float (というかto floatだけすればOKぽい..)
pos = pos.astype(np.float32)
neg = neg.astype(np.float32)

#now we have the large pos and neg file
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.naive_bayes import MultinomialNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier

fold=10
frac_test=0.3



#1. compare tss vs bas only vs bin only
y_real_all = []
y_prob_all = []
y_prob_tss = []
y_prob_bin = []
y_prob_bas = []
for i in range(fold):
    np.random.seed(fold)
    train_pos = np.random.choice(2, pos.shape[0], p=[frac_test, 1-frac_test])
    neg_X_ds = neg.sample(pos.shape[0])  # downsampled so that neg vs pos = 1:1
    train_neg = np.random.choice(2, neg_X_ds.shape[0], p=[frac_test, 1-frac_test])
    X_train = pd.concat([pos[train_pos == 1], neg_X_ds[train_neg == 1]])
    X_test = pd.concat([pos[train_pos == 0], neg_X_ds[train_neg == 0]])
    y_train = np.array([1] * sum(train_pos == 1) + [0] * sum(train_neg == 1))
    y_test = np.array([1] * sum(train_pos == 0) + [0] * sum(train_neg == 0))
    #using all features:
    regr = RandomForestClassifier(random_state=1,n_jobs=-1)
    regr.fit(X_train.fillna(0), y_train)
    y_prob = regr.predict_proba(X_test.fillna(0))
    y_prob = y_prob[:,1]
    y_real_all = y_real_all + list(y_test)
    y_prob_all = y_prob_all + list(y_prob)
    
    #using TSS only (worked):
    regr.fit(np.array(X_train["tss_distance"]).reshape(-1, 1), y_train)
    y_prob = regr.predict_proba(np.array(X_test["tss_distance"]).reshape(-1, 1))
    y_prob = y_prob[:,1]
    y_prob_tss = y_prob_tss + list(y_prob)
    
    #using bin only:
    regr.fit(X_train[list(feat_to_use_bin) + list(feat_to_use_rm)].fillna(0), y_train)
    y_prob = regr.predict_proba(X_test[list(feat_to_use_bin) + list(feat_to_use_rm)].fillna(0))
    y_prob = y_prob[:,1]
    y_prob_bin = y_prob_bin + list(y_prob)
    #using bas only (worked):
    regr.fit(X_train[feat_to_use_bas], y_train)
    y_prob = regr.predict_proba(X_test[feat_to_use_bas])
    y_prob = y_prob[:,1]
    y_prob_bas = y_prob_bas + list(y_prob)

#output as df:
df = pd.DataFrame({"label":y_real_all,"full":y_prob_all, "bin":y_prob_bin, "bas":y_prob_bas, "tss":y_prob_tss})
df.to_csv("rf_evaluation_withinfeat.tsv", sep="\t")



#vs other algorithms
y_real_all = []
y_prob_all_logr = []
y_prob_all_svm = []
y_prob_all_nb = []
y_prob_all_knn = []
y_prob_all_rf = []
y_prob_all_rf_best = []
y_prob_all_adab = []
y_prob_all_gb = []

#also time
times = [] #それぞれ、first iterationだけ

for i in range(fold):
    np.random.seed(fold)
    train_pos = np.random.choice(2, pos.shape[0], p=[frac_test, 1-frac_test])
    neg_X_ds = neg.sample(pos.shape[0])  # downsampled so that neg vs pos = 1:1
    train_neg = np.random.choice(2, neg_X_ds.shape[0], p=[frac_test, 1-frac_test])
    X_train = pd.concat([pos[train_pos == 1], neg_X_ds[train_neg == 1]])
    X_test = pd.concat([pos[train_pos == 0], neg_X_ds[train_neg == 0]])
    y_train = np.array([1] * sum(train_pos == 1) + [0] * sum(train_neg == 1))
    y_test = np.array([1] * sum(train_pos == 0) + [0] * sum(train_neg == 0))
    # train the predictor, rf, best
    regr = RandomForestClassifier(n_estimators=1299, min_samples_split=6,
                                  max_features=0.288, criterion="entropy", random_state=1,n_jobs=-1)
    t1 = tm.time()
    regr.fit(X_train, y_train)
    y_prob = regr.predict_proba(X_test)
    times.append(tm.time()-t1)
    y_prob = y_prob[:,1]
    y_real_all = y_real_all + list(y_test)
    y_prob_all_rf_best = y_prob_all_rf_best + list(y_prob)
    # train the predictor, rf
    regr = RandomForestClassifier(random_state=1,n_jobs=-1)
    t1 = tm.time()
    regr.fit(X_train, y_train)
    y_prob = regr.predict_proba(X_test)
    times.append(tm.time()-t1)
    y_prob = y_prob[:,1]
    y_prob_all_rf = y_prob_all_rf + list(y_prob)
    #adaboost -- 長すぎ...
    regr = AdaBoostClassifier(random_state=1)
    t1 = tm.time()
    regr.fit(X_train, y_train)
    y_prob = regr.predict_proba(X_test)
    times.append(tm.time() - t1)
    y_prob = y_prob[:,1]
    y_prob_all_adab = y_prob_all_adab + list(y_prob)
    #gradient boost -- 多分長い.. けど頑張るしかないぽいな..
    regr = GradientBoostingClassifier(random_state=1)
    t1 = tm.time()
    regr.fit(X_train, y_train)
    y_prob = regr.predict_proba(X_test)
    times.append(tm.time() - t1)
    y_prob = y_prob[:,1]
    y_prob_all_gb = y_prob_all_gb + list(y_prob)
    #logr
    logist_regr = LogisticRegression(random_state=0)
    t1 = tm.time()
    logist_regr.fit(X_train, y_train)
    y_prob = logist_regr.predict_proba(X_test)
    times.append(tm.time() - t1)
    y_prob = y_prob[:, 1]
    y_prob_all_logr = y_prob_all_logr + list(y_prob)
    #svm
    clf = svm.SVC(gamma='auto', probability=True)
    t1 = tm.time()
    clf.fit(X_train, y_train)
    y_prob = clf.predict_proba(X_test)
    times.append(tm.time() - t1)
    y_prob = y_prob[:, 1]
    y_prob_all_svm = y_prob_all_svm + list(y_prob)
    #naive bayes
    nb = MultinomialNB()
    t1 = tm.time()
    nb.fit(X_train, y_train)
    y_prob = nb.predict_proba(X_test)
    times.append(tm.time() - t1)
    y_prob = y_prob[:, 1]
    y_prob_all_nb = y_prob_all_nb + list(y_prob)
    #knn
    knn = KNeighborsClassifier()
    t1 = tm.time()
    knn.fit(X_train, y_train)
    y_prob = knn.predict_proba(X_test)
    times.append(tm.time() - t1)
    y_prob = y_prob[:, 1]
    y_prob_all_knn = y_prob_all_knn + list(y_prob)
    print ("done {0}, {1}".format(i, tm.ctime()))

#output as df:
df = pd.DataFrame({"label":y_real_all,"rf, tuned":y_prob_all_rf_best, "rf":y_prob_all_rf, "adab":y_prob_all_adab,
                   "grad":y_prob_all_gb, "logr":y_prob_all_logr, "svm":y_prob_all_svm, "bayes":y_prob_all_nb,
                   "knn":y_prob_all_knn})
df.to_csv("rf_evaluation_algorithms.tsv", sep="\t")

#and the time
df = pd.DataFrame({"1":times[:7], "2":times[7:7*2], "3":times[7*2:7*3], "4":times[7*3:7*4],"5": times[7*4:7*5], "6":times[7*5:7*6], "7":times[7*6:7*7], "8":times[7*7:7*8], "9":times[7*8:7*9], "10":times[7*9:7*10]})
df.index = ["rf","adab","grad","logr","svm","bayes","knn"]
df.to_csv("rf_evaluation_algorithms_times.tsv", sep="\t")
#and in local:
#no, in script
import subprocess
subprocess.call(["gsutil","cp","rf_evaluation_withinfeat.tsv","gs://qingbowang/"])
subprocess.call(["gsutil","cp","rf_evaluation_algorithms.tsv","gs://qingbowang/"])
subprocess.call(["gsutil","cp","rf_evaluation_algorithms_times.tsv","gs://qingbowang/"])
# gsutil cp rf_evaluation_withinfeat.tsv gs://qingbowang/
# gsutil cp rf_evaluation_algorithms.tsv gs://qingbowang/
# gsutil cp rf_evaluation_algorithms_times.tsv gs://qingbowang/


#その前にこっちでも軽くrocを確認しようか: いやいいや面倒
#もしこれでrfがいい感じじゃなかったらfeature削った後で比べる
fpr, tpr, threshold = metrics.roc_curve(y_real_all, y_prob_all)
roc_auc = metrics.auc(fpr, tpr)

#localかnotebookどっちでもいいので明日やる


