# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

from joblib import dump, load
import pandas as pd
import numpy as np
import sys
import time as tm
import subprocess

tissue_name = sys.argv[1]

print (tm.ctime())
vg = pd.read_csv("vgpair_and_funcannot_{0}.tsv.bgz".format(tissue_name), sep="\t", compression="gzip", engine="c")
#load basenji annotations
dfbas = []
for chk in ["0to499", "500to999", "1000to1499", "1500to1999", "2000to2499", "2500to2999",
            "3000to3499", "3500to3999", "4000to4499", "4500to4999", "5000to5311"]:
    df = pd.read_csv("sad_vartouse_col{0}_{1}.tsv.bgz".format(chk, tissue_name), sep="\t", compression="gzip",
                        engine="c", index_col=0, nrows=10) #just to set the dtype
    float_cols = [c for c in df if df[c].dtype == "float64"]
    float32_cols = {c: np.float32 for c in float_cols}

    df = pd.read_csv("sad_vartouse_col{0}_{1}.tsv.bgz".format(chk, tissue_name), sep="\t", compression="gzip",
                        engine="c", index_col=0, dtype=float32_cols)
    dfbas.append(df)
    print ("done {0}, {1}".format(chk, tm.ctime()))
dfbas = pd.concat(dfbas, axis=1)
print ("before removing duplicates, basenji shape {0}, {1}".format(dfbas.shape, tm.ctime()))
dfbas = dfbas[~dfbas.index.duplicated(keep='first')]
print ("done removing duplicates, basenji shape {0}, {1}".format(dfbas.shape, tm.ctime()))

#load roadmap annotations
dfrm = []
for h in ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3"]:
    df = pd.read_csv("roadmap_vartouse_{0}_{1}.tsv.bgz".format(h, tissue_name), sep="\t", compression="gzip",
                     engine="c", index_col=0)
    #df = pd.read_csv("roadmap_vartouse_{0}_{1}.tsv.bgz".format(h, tissue_name), sep="\t", compression="gzip",
    #                 engine="c", index_col=0, nrows=1000000) #unit test
    dfrm.append(df)
    print("done {0}, {1}".format(h, tm.ctime()))
dfrm = pd.concat(dfrm, axis=1)

vg.index = vg.variant_id.str.replace("_b38","")
dfrm.index = dfrm.index.str.replace("_b38","")
indices = vg.variant_id + "_" + vg.gene_id
feat_to_use_bin = ['DHS_Trynka', 'UTR_3_UCSC', 'UTR_5_UCSC', 'Coding_UCSC',
       'Intron_UCSC', 'TFBS_ENCODE', 'Promoter_UCSC', 'Enhancer_Hoffman'] #all the non-histone feature
vg = vg[feat_to_use_bin + ["tss_distance"]]

print ("starting to join, {0}".format(tm.ctime()))
print (vg.shape)


#and concat them, and also add the basenji features
dfall = vg.join(dfbas, how="left")
print ("done basenji join, {0}".format(tm.ctime()))
print (dfall.shape)

del vg
del dfbas

cname = dfall.columns[-1]
nonnan = (~np.isnan(np.array(dfall[cname])))
indices = indices[nonnan]
dfall = dfall[nonnan]
print ("done removing NAs, {0}".format(tm.ctime()))
print (dfall.shape)

dfall = dfall.join(dfrm, how="left")
print ("done binannot join, {0}".format(tm.ctime()))
print (dfall.shape)
del dfrm


#and turn to log 1p
def abs_log1p(x):
    return (np.log2(abs(x) + 1))

#let the order of the features match
features_ordered_df = pd.read_csv("rf_feat_to_use_ordered.tsv", sep="\t", index_col=0)
features_ordered = features_ordered_df[tissue_name].dropna()
print ("orders changing: {0}".format(sum(dfall.columns != features_ordered)))#to make sure that the order is actually changing
print ("mismatch after sorting: {0}".format(sum(dfall.columns.sort_values() != features_ordered.sort_values())))
dfall = dfall[features_ordered]

#from here, output the scores little by little to save the memory
N = dfall.shape[0]
chksize = 10000000
regr = load("ems_rfmodel_{0}.sav".format(tissue_name))
for i in range(int(np.ceil(N/chksize))):
    print ("starting chunk {0}, {1}".format(i, tm.ctime()))
    df = (dfall.iloc[chksize*i:chksize*(i+1),:]*1.0).astype(float).apply(lambda x: abs_log1p(x)) #turning to log

    print ("done taking abs, {0}".format(tm.ctime()))

    #and score
    y_prob = regr.predict_proba(df)

    print ("done scoring, {0}".format(tm.ctime()))

    y_prob = y_prob[:,1]
    y_prob = pd.Series(y_prob)
    y_prob.index = indices[chksize*i:chksize*(i+1)]
    y_prob.to_csv("ems_rawscore_gtexvg_all_{0}_chunk{1}.tsv.gz".format(tissue_name, i), sep="\t", compression="gzip")

    print ("done to csv, {0}".format(tm.ctime()))

    #and cp it to google cloud
    import subprocess
    subprocess.call(["gsutil", "cp" ,"./ems_rawscore_gtexvg_all_{0}_chunk{1}.tsv.gz".format(tissue_name, i) ,"gs://qingbowang/ems_v1_test/ems_rawscore_gtexvg_all{0}_chunk{1}.tsv.gz".format(tissue_name, i)])

#finally, remove all the files from local
subprocess.call(["rm", "*tsv.gz"])
subprocess.call(["rm", "*tsv.bgz"])
subprocess.call(["rm", "*.sav"])
