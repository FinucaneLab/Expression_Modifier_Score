# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

import pandas as pd
import sys
import time as tm
import numpy as np
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')



#format the annotated GTEx files and export as tsv.gz

bas_feat = hl.import_table("gs://qingbowang/ems_v1_test/bas_feature_top100_touse.tsv", impute=True).to_pandas()
bin_feat = hl.import_table("gs://qingbowang/ems_v1_test/roadmapannot_feat_to_use_all.tsv", impute=True).to_pandas()

feat_to_use_bin = ['DHS_Trynka', 'UTR_3_UCSC', 'UTR_5_UCSC', 'Coding_UCSC',
       'Intron_UCSC', 'TFBS_ENCODE', 'Promoter_UCSC', 'Enhancer_Hoffman'] #all the non-histone features

fs = ["gs://gnomad-qingbowang/sad_allvar_hg38_col0to499_wID.ht",
"gs://gnomad-qingbowang/sad_allvar_hg38_col500to999_wID.ht",
"gs://gnomad-qingbowang/sad_allvar_hg38_col1000to1499_wID.ht",
"gs://gnomad-qingbowang/sad_allvar_hg38_col1500to1999_wID.ht",
"gs://gnomad-qingbowang/sad_allvar_hg38_col2000to2499_wID.ht",
"gs://gnomad-qingbowang/sad_allvar_hg38_col2500to2999_wID.ht",
"gs://gnomad-qingbowang/sad_allvar_hg38_col3000to3499_wID.ht",
"gs://gnomad-qingbowang/sad_allvar_hg38_col3500to3999_wID.ht",
"gs://gnomad-qingbowang/sad_allvar_hg38_col4000to4499_wID.ht",
"gs://gnomad-qingbowang/sad_allvar_hg38_col4500to4999_wID.ht",
"gs://gnomad-qingbowang/sad_allvar_hg38_col5000to5311_wID.ht"]
fs_num = pd.Series(fs).apply(lambda x: x.split("_")[-2])


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



for tissue_name in tissues:
    print("start {0}, {1}".format(tissue_name, tm.ctime()))

    feat_to_use_bas = list(bas_feat[tissue_name])
    feat_to_use_roadmap = list(bin_feat[tissue_name].dropna())
    #export basenji features
    i = 0 #for cnt
    for f in fs:
        htsub = hl.read_table(f).repartition(100) #smaller partition solves everything? seems like it does make things better
        feat_to_use_mock_sub = np.intersect1d(list(htsub.row), feat_to_use_bas)
        htsub = htsub.select(*feat_to_use_mock_sub)
        htsub.export("gs://qingbowang/ems_v1_test/bgz_chunks/sad_vartouse_{0}_{1}.tsv.bgz".format(fs_num[i], tissue_name))
        i += 1
    #tissue specific bin annots:
    marks = ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3" ,"H3K9ac", "H3K9me3"]
    for h in marks:
        htsub = hl.read_table("gs://qingbowang/gtex_v8_wgs_binannot_pertissue_{0}.ht".format(h))
        feat_to_use_mock_sub = np.intersect1d(list(htsub.row), feat_to_use_roadmap)
        htsub = htsub.select(*feat_to_use_mock_sub)
        htsub.export("gs://qingbowang/ems_v1_test/bgz_chunks/roadmap_vartouse_{0}_{1}.tsv.bgz".format(h, tissue_name))

    print ("done {0}, {1}".format(tissue_name, tm.ctime()))

    #also export the binary annotations, and tss distance
    vg = hl.read_table("gs://qingbowang/ems_v1_test/{0}_allpairs_binannot_fmannot.ht".format(tissue_name))
    feats = ["tss_distance"] + feat_to_use_bin #variant and gene id are keys
    vg.select(*feats).export("gs://qingbowang/ems_v1_test/bgz_chunks/vgpair_and_funcannot_{0}.tsv.bgz".format(tissue_name))
