# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

import sys
import pandas as pd
import numpy as np
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')
#also changed this to be done in hail cluster rather than single vm

#collect all the bas features, and take top 100
chr_list = []
for i in range(1,22+1):
    chr_list.append("chr" + str(i))
chr_list.append("chrX")
tissue_name = "Whole_Blood"


imps = []
feat_to_use = {}
for chr in chr_list: #chr to leave out
    fn = "gs://qingbowang/ems_v1_test/bas_feat_imp_100trees_{0}_{1}_leftout.tsv".format(tissue_name, chr)
    with hl.hadoop_open(fn, 'r') as f:
        imp = pd.read_csv(f, sep="\t", index_col=0, header=None)#これがdownstreamを全て狂わせた.. hadoopに移ったのでheader=None必要..
    imps.append(imp) #the full list of importances
    feat_to_use[chr] = list(imp.iloc[:100,:].index) #top 100 are the features to use

imps = pd.concat(imps, axis=1, sort=True) #check whether the axis is right / the order actually changes
feat_to_use = pd.DataFrame(feat_to_use)
feat_to_use.columns = chr_list

fn = "gs://qingbowang/ems_v1_test/bas_feature_top100_touse_per_leftout_chr.tsv"
with hl.hadoop_open(fn, 'w') as f:
    feat_to_use.to_csv(f, sep="\t")
fn = "gs://qingbowang/ems_v1_test/bas_feature_imp_100forest_per_leftout_chr.tsv"
with hl.hadoop_open(fn, 'w') as f:
    imps.to_csv(f, sep="\t")
