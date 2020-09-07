# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

#collect all the bas features, and take top 100
tissue_list = ["Whole_Blood",
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
import pandas as pd
imps = []
feat_to_use = {}
for tissue_name in tissue_list:
    imp = pd.read_csv("bas_feature_imp_100forest_{0}.tsv".format(tissue_name), sep="\t", index_col=0)
    imps.append(imp) #the full list of importances
    feat_to_use[tissue_name] = list(imp.iloc[:100,:].index) #top 100 are the features to use

imps = pd.concat(imps, axis=1) #check whether the axis is right / the order actually changes
feat_to_use = pd.DataFrame(feat_to_use)
feat_to_use.columns = tissue_list
feat_to_use.to_csv("bas_feature_top100_touse.tsv", sep="\t")
imps.to_csv("bas_feature_imp_100forest_alltissues.tsv", sep="\t")




