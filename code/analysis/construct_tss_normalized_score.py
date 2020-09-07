# -*- coding: utf-8 -*-
__author__ = 'QingboWang'
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')
import pandas as pd
import time as tm


tissues = ["Whole_Blood",
            "Muscle_Skeletal",
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
            "Liver",
            "Adipose_Visceral_Omentum",
            "Pancreas",
            "Lung",
            "Pituitary",
            "Brain_Nucleus_accumbens_basal_ganglia",
            "Colon_Transverse",
            "Adipose_Subcutaneous",
            "Esophagus_Mucosa",
            "Brain_Cerebellum",
            "Brain_Cortex",
            "Thyroid",
            "Stomach",
            "Breast_Mammary_Tissue",
            "Colon_Sigmoid",
            "Skin_Not_Sun_Exposed_Suprapubic",
            "Testis"]

for tissue_name in tissues:
    print ("starting {0}, {1}".format(tissue_name, tm.ctime()))
    ht = hl.read_table("gs://qingbowang/ems_v1_test/ems_pcausal_gtexvg_all{0}.ht".format(tissue_name))
    ht = ht.key_by("vg")
    #annotate the TSS distance
    vg = hl.read_table("gs://qingbowang/ems_v1_test/{0}_allpairs.ht".format(tissue_name))
    vg = vg.annotate(vg=vg.variant_id + "_" + vg.gene_id)
    vg = vg.key_by("vg")
    ht = ht.annotate(abs_tss_dist = hl.abs(vg[ht.key].tss_distance))
    agged = ht.group_by(ht.abs_tss_dist).aggregate(n=hl.agg.count(), ave=hl.agg.mean(ht.p_causal))
    agged.write("gs://qingbowang/ems_v1_test/{0}_per_tss_dist_p_causal.ht".format(tissue_name), overwrite=True)

#just for whole blood, smooth it
tissue_name = "Whole_Blood"
ht = hl.read_table("gs://qingbowang/ems_v1_test/{0}_per_tss_dist_p_causal.ht".format(tissue_name))
df = ht.to_pandas()
df["weighted_sum"] = df.n * df.ave
n_smoothed = df.n.rolling(100, min_periods=1).sum()
weighted_sum_smoothed = df.weighted_sum.rolling(100, min_periods=1).sum()
df["smoothed_ave"] = weighted_sum_smoothed / n_smoothed
hl.Table.from_pandas(df[["abs_tss_dist", "smoothed_ave"]]).write("gs://qingbowang/ems_v1_test/{0}_per_tss_dist_p_causal_smoothed.ht".format(tissue_name))




