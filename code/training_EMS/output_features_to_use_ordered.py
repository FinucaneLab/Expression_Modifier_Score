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
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/n_pos_and_neg_per_tissue.tsv", "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/bas_feature_top100_touse.tsv", "./"])
subprocess.call(["gsutil", "cp" ,"gs://qingbowang/ems_v1_test/roadmapannot_feat_to_use_all.tsv", "./"])
basfeats = pd.read_csv("bas_feature_top100_touse.tsv", sep="\t", index_col=0)
rmfeats = pd.read_csv("roadmapannot_feat_to_use_all.tsv", sep="\t", index_col=0)

#load the data (load from the cloud first)
import time as tm

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

feat_to_use = {}
for tissue_name in tissues:

    print ("starting {0}, {1}".format(tissue_name, tm.ctime()))

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



    #subset to features to use
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
    #and turn the log1p この文脈ではno needだけど一応
    def abs_log1p(x):
        return (np.log2(abs(x) + 1))
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

    feat_to_use[tissue_name] = list(pos.columns)

    #and remove things in local
    subprocess.call(["rm" ,"grid_search_result_{0}.tsv".format(tissue_name)])
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



pd.DataFrame.from_dict(feat_to_use, orient="index").T.to_csv("rf_feat_to_use_ordered.tsv", sep="\t")
subprocess.call(["gsutil", "cp" ,"rf_feat_to_use_ordered.tsv", "gs://qingbowang/ems_v1_test/rf_feat_to_use_ordered.tsv"])