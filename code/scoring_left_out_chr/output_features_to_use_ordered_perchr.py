# -*- coding: utf-8 -*-
__author__ = 'QingboWang'
import sys
import pandas as pd
import numpy as np
import time as tm
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import subprocess


#non tissue specific things
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/n_pos_and_neg_per_leftout_chr.tsv", "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/bas_feature_top100_touse_per_leftout_chr.tsv", "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/roadmapannot_feat_to_use_per_leftout_chr.tsv", "./"])
basfeats = pd.read_csv("bas_feature_top100_touse_per_leftout_chr.tsv", sep="\t", index_col=0)
rmfeats = pd.read_csv("roadmapannot_feat_to_use_per_leftout_chr.tsv", sep="\t", index_col=0)

#load the data (load from the cloud first)
import time as tm
tissue_name = "Whole_Blood"
chr_list = []
for i in range(1,22+1):
    chr_list.append("chr" + str(i))
chr_list.append("chrX")

feat_to_use = {}

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

posbas0 = pd.read_csv("{0}_positive_training_vg_annotated.tsv".format(tissue_name), sep="\t", index_col=0)
negbas0 = pd.read_csv("{0}_negative_training_vg_annotated.tsv".format(tissue_name), sep="\t", index_col=0)

marks = ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3"]
posrm0 = {}
negrm0 = {}
for m in marks:
    posrm0[m] = pd.read_csv("{0}_positive_training_vg_roadmapannot_{1}.tsv".format(tissue_name, m), sep="\t", index_col=0)
    negrm0[m] = pd.read_csv("{0}_negative_training_vg_roadmapannot_{1}.tsv".format(tissue_name, m), sep="\t", index_col=0)

for chr in chr_list:

    print ("starting {0}, {1}".format(chr, tm.ctime()))

    #subset to features to use
    feat_to_use_bas = basfeats.loc[:,chr]
    feat_to_use_rm = rmfeats.loc[:,chr].dropna()
    feat_to_use_bin = ['DHS_Trynka', 'UTR_3_UCSC', 'UTR_5_UCSC', 'Coding_UCSC',
           'Intron_UCSC', 'TFBS_ENCODE', 'Promoter_UCSC', 'Enhancer_Hoffman'] #all the non-histone features

    #load data
    pos = []
    neg = []
    #filter out the left out chr
    posbas = posbas0[posbas0.index.str.split("_").str[0] != chr]
    negbas = negbas0[negbas0.index.str.split("_").str[0] != chr]

    posbas = posbas.loc[:,feat_to_use_bin + list(feat_to_use_bas)]
    negbas = negbas.loc[:,feat_to_use_bin + list(feat_to_use_bas)]
    def abs_log1p(x):
        return (np.log2(abs(x) + 1))
    posbas = posbas.apply(lambda x: abs_log1p(x))
    negbas = negbas.apply(lambda x: abs_log1p(x))

    pos.append(posbas)
    neg.append(negbas)
    marks = ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3" ,"H3K9ac", "H3K9me3"]
    for m in marks:
        # filter to leave one chr
        posrm = posrm0[m][posrm0[m].index.str.split("_").str[0] != chr]
        negrm = negrm0[m][negrm0[m].index.str.split("_").str[0] != chr]

        posrm = posrm.loc[:, feat_to_use_rm[feat_to_use_rm.str.contains(m)]]
        negrm = negrm.loc[:, feat_to_use_rm[feat_to_use_rm.str.contains(m)]]
        posrm = posrm*1
        negrm = negrm*1
        pos.append(posrm)
        neg.append(negrm)
    pos = pd.concat(pos, axis=1)
    neg = pd.concat(neg, axis=1)

    feat_to_use[chr] = list(pos.columns)


pd.DataFrame.from_dict(feat_to_use, orient="index").T.to_csv("rf_feat_to_use_ordered_per_leftout_chr.tsv", sep="\t")
subprocess.call(["gsutil", "cp" ,"rf_feat_to_use_ordered_per_leftout_chr.tsv", "gs://qingbowang/ems_v1_test/rf_feat_to_use_ordered_per_leftout_chr.tsv"])


#and remove things in local
subprocess.call(["rm" ,"{0}_positive_training_vg_annotated.tsv".format(tissue_name)])
subprocess.call(["rm" ,"{0}_negative_training_vg_annotated.tsv".format(tissue_name)])

subprocess.call(["rm" ,"{0}_negative_training_vg_roadmapannot_H3K27ac.tsv".format(tissue_name)])
subprocess.call(["rm" ,"{0}_negative_training_vg_roadmapannot_H3K27me3.tsv".format(tissue_name)])
subprocess.call(["rm" ,"{0}_negative_training_vg_roadmapannot_H3K36me3.tsv".format(tissue_name)])
subprocess.call(["rm" ,"{0}_negative_training_vg_roadmapannot_H3K4me1.tsv".format(tissue_name)])
subprocess.call(["rm" ,"{0}_negative_training_vg_roadmapannot_H3K4me3.tsv".format(tissue_name)])
subprocess.call(["rm" ,"{0}_negative_training_vg_roadmapannot_H3K9ac.tsv".format(tissue_name)])
subprocess.call(["rm" ,"{0}_negative_training_vg_roadmapannot_H3K9me3.tsv".format(tissue_name)])

subprocess.call(["rm" ,"{0}_positive_training_vg_roadmapannot_H3K27ac.tsv".format(tissue_name)])
subprocess.call(["rm" ,"{0}_positive_training_vg_roadmapannot_H3K27me3.tsv".format(tissue_name)])
subprocess.call(["rm" ,"{0}_positive_training_vg_roadmapannot_H3K36me3.tsv".format(tissue_name)])
subprocess.call(["rm" ,"{0}_positive_training_vg_roadmapannot_H3K4me1.tsv".format(tissue_name)])
subprocess.call(["rm" ,"{0}_positive_training_vg_roadmapannot_H3K4me3.tsv".format(tissue_name)])
subprocess.call(["rm" ,"{0}_positive_training_vg_roadmapannot_H3K9ac.tsv".format(tissue_name)])
subprocess.call(["rm" ,"{0}_positive_training_vg_roadmapannot_H3K9me3.tsv".format(tissue_name)])



