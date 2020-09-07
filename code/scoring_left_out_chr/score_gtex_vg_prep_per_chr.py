# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

import pandas as pd
import sys
import time as tm
import numpy as np
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')

bas_feat = hl.import_table("gs://qingbowang/ems_v1_test/bas_feature_top100_touse_per_leftout_chr.tsv", impute=True).to_pandas()
bin_feat = hl.import_table("gs://qingbowang/ems_v1_test/roadmapannot_feat_to_use_per_leftout_chr.tsv", impute=True).to_pandas()

feat_to_use_bin = ['DHS_Trynka', 'UTR_3_UCSC', 'UTR_5_UCSC', 'Coding_UCSC',
       'Intron_UCSC', 'TFBS_ENCODE', 'Promoter_UCSC', 'Enhancer_Hoffman'] #all the non-histone features

fs = ["gs://qingbowang/basenji_scores/sad_allvar_hg38_col0to499_wID.ht",
"gs://qingbowang/basenji_scores/sad_allvar_hg38_col500to999_wID.ht",
"gs://qingbowang/basenji_scores/sad_allvar_hg38_col1000to1499_wID.ht",
"gs://qingbowang/basenji_scores/sad_allvar_hg38_col1500to1999_wID.ht",
"gs://qingbowang/basenji_scores/sad_allvar_hg38_col2000to2499_wID.ht",
"gs://qingbowang/basenji_scores/sad_allvar_hg38_col2500to2999_wID.ht",
"gs://qingbowang/basenji_scores/sad_allvar_hg38_col3000to3499_wID.ht",
"gs://qingbowang/basenji_scores/sad_allvar_hg38_col3500to3999_wID.ht",
"gs://qingbowang/basenji_scores/sad_allvar_hg38_col4000to4499_wID.ht",
"gs://qingbowang/basenji_scores/sad_allvar_hg38_col4500to4999_wID.ht",
"gs://qingbowang/basenji_scores/sad_allvar_hg38_col5000to5311_wID.ht"] #ugly but works...
fs_num = pd.Series(fs).apply(lambda x: x.split("_")[-2])

tissue_name = "Whole_Blood"
chr_list = []
for i in range(1,22+1):
    chr_list.append("chr" + str(i))
chr_list.append("chrX")

for chr in chr_list:
    print("start {0}, {1}".format(chr, tm.ctime()))

    feat_to_use_bas = list(bas_feat[chr])
    feat_to_use_roadmap = list(bin_feat[chr].dropna())
    i = 0 #for cnt
    for f in fs:
        htsub = hl.read_table(f)#.repartition(100) #smaller partition solves everything? seems like it does make things better

        #filter to the chr of interest
        htsub = htsub.filter(htsub.ID.split("_")[0]==chr)

        feat_to_use_mock_sub = np.intersect1d(list(htsub.row), feat_to_use_bas)
        htsub = htsub.select(*feat_to_use_mock_sub)
        htsub.export("gs://qingbowang/ems_v1_test/bgz_chunks/sad_vartouse_{0}_{1}.tsv.bgz".format(fs_num[i], chr))
        i += 1
    #tissue specific bin annots:
    marks = ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3" ,"H3K9ac", "H3K9me3"]
    for h in marks:
        htsub = hl.read_table("gs://qingbowang/gtex_v8_wgs_binannot_pertissue_{0}.ht".format(h))

        #filter to the chr of interest
        htsub = htsub.filter(htsub.variant_id.split("_")[0]==chr)

        feat_to_use_mock_sub = np.intersect1d(list(htsub.row), feat_to_use_roadmap)
        htsub = htsub.select(*feat_to_use_mock_sub)
        htsub.export("gs://qingbowang/ems_v1_test/bgz_chunks/roadmap_vartouse_{0}_{1}.tsv.bgz".format(h, chr))

    print ("done {0}, {1}".format(chr, tm.ctime()))


    #also export the binary annotations, and tss distance, per chr as well
    vg = hl.read_table("gs://qingbowang/ems_v1_test/{0}_allpairs_binannot_fmannot.ht".format(tissue_name))#bin annotはここでしてあるからokay
    feats = ["tss_distance"] + feat_to_use_bin #variant and gene id are keys
    vg.filter(vg.variant_id.split("_")[0]==chr).select(*feats).export("gs://qingbowang/ems_v1_test/bgz_chunks/vgpair_and_funcannot_{0}.tsv.bgz".format(chr))


#and export everything in a single huge VM, do everything below in the VM
