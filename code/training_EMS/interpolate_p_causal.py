# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

#do this in random vm
import subprocess
import pandas as pd
import numpy as np

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

r = pd.DataFrame(pd.Series(np.arange(0,1.000001, 0.001)))
r.columns = ["rf_score_bin"]
r.index = r.rf_score_bin
for tissue_name in tissues:
    fn = "gs://qingbowang/ems_v1_test/ems_p_causal_{0}.tsv".format(tissue_name)
    subprocess.call(["gsutil", "cp", fn, "./"])
    pcausal = pd.read_csv("ems_p_causal_{0}.tsv".format(tissue_name), sep="\t", index_col=0)
    j = r.join(pcausal).interpolate(limit_direction="both") #joined, and interpolated for nans
    j["confidence_gain_log10"] = np.log10(j.confidence_gain)#need to re-calculate the log10 for interpolated ones
    #check there's no NA
    print("num. NAs: {0}".format(j.isna().sum().sum()))
    #also sanity check
    if tissue_name=="Whole_Blood":
        print (j.head())
        print (j.tail())

    if tissue_name=="Heart_Atrial_Appendage":
        print (j.head())
        print (j.tail())
    j.to_csv("ems_p_causal_interpolated_{0}.tsv".format(tissue_name), sep="\t")
    subprocess.call(["gsutil", "cp", "ems_p_causal_interpolated_{0}.tsv".format(tissue_name), "gs://qingbowang/ems_v1_test/"])
