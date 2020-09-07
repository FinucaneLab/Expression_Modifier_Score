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

import pickle
from joblib import dump, load
tissue_name = "Whole_Blood"
regr = load("/Users/qingbowang/Downloads/ems_rfmodel_{0}.sav".format(tissue_name))
features_ordered_df = pd.read_csv("~/Downloads/rf_feat_to_use_ordered.tsv", sep="\t", index_col=0)

#btw, num. features:
n_feat = features_ordered_df.apply(lambda x: sum(~x.isnull()), axis=0)

regr = load("./ems_rfmodel_{0}.sav".format(tissue_name))
features_ordered_df = pd.read_csv("./rf_feat_to_use_ordered.tsv", sep="\t", index_col=0)

importances = regr.feature_importances_

imp = pd.DataFrame(importances)
imp.columns = ["imp"]
imp.index = list(features_ordered_df.Whole_Blood.dropna())
imp.sort_values(by="imp", ascending=False).to_csv("wb_feat_imp.tsv", sep="\t")

imp = pd.read_csv("~/Downloads/wb_feat_imp.tsv", sep="\t", index_col=0)

#skimming the results manually:
"""
tss_distance	0.4311639186204463

CNhs13058	0.0324490177270843
ENCFF891PZS	0.0296225527736572
CNhs11675	0.01873899011772308
ENCFF628GHA	0.017693406063407
ENCFF634ZUJ	0.012923342725049166

their name: got from supp file:
CAGE:acute myeloid leukemia (FAB M5) cell line
CHIP:H3K4me3:neutrophil male
CAGE:Whole blood (ribopure), , donation1
CHIP:H3K4me3:neutrophil
DNASE:HL-60

SpleenE113-H3K27ac	0.002212202680979495
CD15_Primary_CellsE030-H3K4me3	0.00167114702013523
Monocytes-CD14+_RO01746E124-H3K4me3	0.0011357954076057616
Monocytes-CD14+_RO01746E124-H3K27ac	0.0006320041029799046
CD15_Primary_CellsE030-H3K4me1	0.0005951577180239021

DHS_Trynka	0.004458935389005094
Coding_UCSC	0.0016388569771689824
UTR_3_UCSC	0.0008694423186085235
Intron_UCSC	0.0004260769616154649
TFBS_ENCODE	0.0004127870367903389
"""

#also want the total %:
imp = pd.read_csv("~/Downloads/wb_feat_imp_manannot.tsv", sep="\t", index_col=0)
imp.groupby("type").imp.sum() #this is the per feature type
#prints:
"""
bas                 0.549883
tisssue_spec_bin    0.010224
tss                 0.431164
union_bin           0.008729
"""

#also annotate their real name as well, manually
df = pd.read_csv("~/Downloads/feat_imp_for_fig2.txt", sep="\t")
df["cols"] = ["tab:purple"]*1 + ["#ffffffff"]*1 + ["tab:green"]*5 + ["#ffffffff"]*2 + ["tab:brown"]*5 + ["#ffffffff"]*2 + ["tab:blue"]*5
x = df.index
y = df.imp * 100
n = df["name"].str.replace("dm5","  ").str.replace("dm4","   ")\
              .str.replace("dm3","    ").str.replace("dm2","     ").str.replace("dm","      ").replace("tss_distance","TSS distance")

import matplotlib.lines as mlines
mk0 = mlines.Line2D([], [], color="tab:purple", marker='s', linestyle='None',
                          markersize=10, label='distance to TSS')
mk1 = mlines.Line2D([], [], color="tab:green", marker='s', linestyle='None',
                          markersize=10, label='Basenji scores')
mk2 = mlines.Line2D([], [], color="tab:brown", marker='s', linestyle='None',
                          markersize=10, label='baseline annotations')
mk3 = mlines.Line2D([], [], color="tab:blue", marker='s', linestyle='None',
                          markersize=10, label='histone marks')

fig, (ax1,ax2) = plt.subplots(2,1,sharex=True,
                         figsize=(8,6))
ax1.spines['bottom'].set_visible(False)
ax1.tick_params(axis='x',which='both',bottom=False)
ax2.spines['top'].set_visible(False)

ax2.set_ylim(0,3.5)
ax1.set_ylim(43.5-3.5,43.5)
bars1 = ax1.bar(n, y, color=list(df.cols))
bars2 = ax2.bar(n, y, color=list(df.cols))

for tick in ax2.get_xticklabels():
    tick.set_rotation(90)
d = .005
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
kwargs.update(transform=ax2.transAxes)
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
ax2.text(7-0.5,0.1,"...", fontsize=18)
ax2.text(14-0.5,0.1,"...", fontsize=18)
ax2.text(21-0.5,0.1,"...", fontsize=18)
plt.tight_layout()
fig.text(0, 0.7, "Feature importance (%)", va='center', rotation='vertical', fontsize=14) #ylabel
plt.xlabel("Features", fontsize=16)
ax2.set_xticks(df[~df["name"].str.contains("dm")].index)
ax1.legend(handles=[mk0,mk1,mk2,mk3], title="Feature category")
plt.savefig("/Users/qingbowang/downloads/fig2_feat_imp.png", dpi=400)


#also the bar chart:
y1 = 0.010224*100
y2 = 0.008729*100
y3 = 0.549883*100
y4 = 0.431164*100
fig = plt.figure(figsize=(5.2,5.2))
plt.bar([1,2,3], [y1,0,0], label='y1', color="tab:blue")
plt.bar([1,2,3], [y2,0,0] ,bottom=y1,label='y2', color="tab:brown")
plt.bar([1,2,3], [y3,0,0] ,bottom=y1+y2,label='y3', color="tab:green")
plt.bar([1,2,3], [y4,0,0] ,bottom=y1+y2+y3,label='y4', color="tab:purple")
plt.ylabel("Feature importance (%) per category", fontsize=16)
import matplotlib.lines as mlines
mk0 = mlines.Line2D([], [], color="tab:purple", marker='s', linestyle='None',
                          markersize=10, label='distance to TSS\n(43.1%, n=1)')
mk1 = mlines.Line2D([], [], color="tab:green", marker='s', linestyle='None',
                          markersize=10, label='Basenji scores\n(55.0%, n=5,313)')
mk2 = mlines.Line2D([], [], color="tab:brown", marker='s', linestyle='None',
                          markersize=10, label='baseline annotations\n(0.87%, n=12)')
mk3 = mlines.Line2D([], [], color="tab:blue", marker='s', linestyle='None',
                          markersize=10, label='histone marks\n(1.02%, n=795)')
plt.legend(handles=[mk0,mk1,mk2,mk3], loc="center right",
           title="Feature category\n(total feature importance,\n and the number of features)")
plt.tick_params(
    axis='x',          
    which='both',      
    bottom=False,      
    top=False,         
    labelbottom=False) 
plt.tight_layout()
plt.savefig("/Users/qingbowang/downloads/fig2_feat_imp_bar.png", dpi=400)


#also this for all tisssues, for the supp file:
#this is done in vm:
tissues = ["Whole_Blood",
            "Muscle_Skeletal",
           "Liver",
           "Brain_Cerebellum",
            "Prostate",
            "Spleen",
            "Skin_Sun_Exposed_Lower_leg",
            "Artery_Coronary",
            "Esophagus_Muscularis",
            "Esophagus_Gastroesophageal_Junction",
            "Artery_Tibial",
            "Heart_Atrial_Appendage",
            "Nerve_Tibial",
            "Heart_Left_Ventricle",
            "Adrenal_Gland",
            "Adipose_Visceral_Omentum",
            "Pancreas",
            "Lung",
            "Pituitary",
            "Brain_Nucleus_accumbens_basal_ganglia",
            "Colon_Transverse",
            "Adipose_Subcutaneous",
            "Esophagus_Mucosa",
            "Brain_Cortex",
            "Thyroid",
            "Stomach",
            "Breast_Mammary_Tissue",
            "Colon_Sigmoid",
            "Skin_Not_Sun_Exposed_Suprapubic",
            "Testis",
            "Artery_Aorta",
            "Brain_Amygdala",
            "Brain_Anterior_cingulate_cortex_BA24",
            "Brain_Caudate_basal_ganglia",
            "Brain_Cerebellar_Hemisphere",
            "Brain_Frontal_Cortex_BA9",
            "Brain_Hippocampus",
            "Brain_Hypothalamus",
            "Brain_Putamen_basal_ganglia",
            "Brain_Spinal_cord_cervical_c-1",
            "Brain_Substantia_nigra",
            "Cells_Cultured_fibroblasts",
            "Cells_EBV-transformed_lymphocytes",
            "Kidney_Cortex",
            "Minor_Salivary_Gland",
            "Ovary",
            "Small_Intestine_Terminal_Ileum",
            "Uterus",
            "Vagina"]

features_ordered_df = pd.read_csv("./rf_feat_to_use_ordered.tsv", sep="\t", index_col=0)
for tissue_name in tissues:
    subprocess.call(["gsutil", "cp", "gs://qingbowang/ems_v1_test/ems_rfmodel_{0}.sav".format(tissue_name), "./"])
    regr = load("./ems_rfmodel_{0}.sav".format(tissue_name))
    importances = regr.feature_importances_
    imp = pd.DataFrame(importances)
    imp.columns = ["importance"]
    imp.index = list(features_ordered_df[tissue_name].dropna())
    imp["feature"] = imp.index
    imp[["feature","importance"]].sort_values(by="importance", ascending=False).to_csv("{0}_feat_imp.tsv".format(tissue_name), sep="\t")
    subprocess.call(["gsutil", "cp","{0}_feat_imp.tsv".format(tissue_name), "gs://qingbowang/tmp/"])
    print ("copied {0}".format(tissue_name))


