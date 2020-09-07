# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')
import pandas as pd
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

for tissue_name in tissues:

    print ("starting {0}, {1}".format(tissue_name, tm.ctime()))

    fn = "gs://qingbowang/ems_v1_test/ems_p_causal_interpolated_{0}.tsv".format(tissue_name)
    with hl.hadoop_open(fn, 'r') as f:
        pcausal = pd.read_csv(f, sep="\t", index_col=0)
    pcausal["rf_score_bin"] = pcausal.index
    del pcausal["rf_score_bin.1"] #duplicated columns
    
    pcausal = hl.Table.from_pandas(pcausal)
    pcausal = pcausal.transmute(rf_score_bin = hl.format('%.3f', pcausal["rf_score_bin"]))
    
    #score all chunks:
    #get the max
    for i in range(10000): #just take the upperbound
        if not hl.hadoop_exists("gs://qingbowang/ems_v1_test/ems_rawscore_gtexvg_all{0}_chunk{1}.tsv.gz".format(tissue_name, i)):
            imax = i
            break
    dfall = []
    for i in range(imax):
        print ("starting chunk {0} of {1}, {2}".format(i, imax-1, tm.ctime()))
        df = hl.import_table("gs://qingbowang/ems_v1_test/ems_rawscore_gtexvg_all{0}_chunk{1}.tsv.gz".format(tissue_name, i), force=True, impute=True)
        df = df.repartition(80) #80 partition for 1000mann lines
        df = df.annotate(rf_score_bin = hl.format('%.3f', df["0"]))
        pcausal = pcausal.key_by("rf_score_bin")
        df = df.key_by("rf_score_bin")
        df = df.join(pcausal, how="left")
        df = df.rename({'' : 'vg', '0' : 'rf_score_raw'})
        df.write("gs://qingbowang/ems_v1_test/tmp/ems_pcausal_gtexvg_all{0}_chunk{1}.ht".format(tissue_name, i), overwrite=True)
        #read for the final concat
        df = hl.read_table("gs://qingbowang/ems_v1_test/tmp/ems_pcausal_gtexvg_all{0}_chunk{1}.ht".format(tissue_name, i))
        dfall.append(df)
    #and finally, combine them all
    dfall = dfall[0].union(*dfall[1:])
    dfall.write("gs://qingbowang/ems_v1_test/ems_pcausal_gtexvg_all{0}.ht".format(tissue_name), overwrite=True)