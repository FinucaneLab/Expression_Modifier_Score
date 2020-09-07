# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

import sys
import pandas as pd
import numpy as np
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')


tissue_name = sys.argv[1]

#3. take the positive and negative training representative, and export for feature selection and training
vg = hl.read_table("gs://qingbowang/ems_v1_test/{0}_allpairs_binannot_fmannot.ht".format(tissue_name))
N_all = vg.count()
pos = vg.filter( (vg.pp_susie>0.9) & (vg.pp_fm>0.9))
n_pos = pos.count()
neg = vg.filter((vg.pp_susie<0.0001) & (vg.pp_fm<0.0001))
n_neg_tmp = neg.count()
neg = neg.sample( min(1, 10000./n_neg_tmp), seed=1) #sampling 10000, as representative (in reality: want all)
#export them as pandas
pos.export("gs://qingbowang/ems_v1_test/{0}_positive_training_vg_beforeannot.tsv".format(tissue_name))
neg.export("gs://qingbowang/ems_v1_test/{0}_negative_training_vg_beforeannot.tsv".format(tissue_name))



#for those pos and neg training sites, also extract all the basenji scores, (in local)
pos = hl.import_table("gs://qingbowang/ems_v1_test/{0}_positive_training_vg_beforeannot.tsv".format(tissue_name))
neg = hl.import_table("gs://qingbowang/ems_v1_test/{0}_negative_training_vg_beforeannot.tsv".format(tissue_name))
pos = pos.annotate(ID = pos.variant_id.replace("_b38","")).key_by("ID") #since we cut this b_38 in the basenji file..
neg = neg.annotate(ID = neg.variant_id.replace("_b38","")).key_by("ID") #since we cut this b_38 in the basenji file..
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
"gs://qingbowang/basenji_scores/sad_allvar_hg38_col5000to5311_wID.ht"] #hail table with 5000 columns crashes, so dividing it to several chunks

i = 0
for f in fs:
    htsub = hl.read_table(f)
    bas_pos = htsub.filter(hl.is_defined(pos[htsub.ID]))
    bas_neg = htsub.filter(hl.is_defined(neg[htsub.ID]))
    bas_pos.export("gs://qingbowang/ems_v1_test/{0}_sad_postraining_colchunk{1}.tsv".format(tissue_name, i))
    bas_neg.export("gs://qingbowang/ems_v1_test/{0}_sad_negtraining_colchunk{1}.tsv".format(tissue_name, i))
    i += 1

#4. collect those pandas pieces, and do the feature selection and training

#4.1. collect those pieces, annotate v-g with basenji scores by join
cols_to_use_bin = ['DHS_Trynka', 'UTR_3_UCSC', 'UTR_5_UCSC', 'Coding_UCSC', 'Intron_UCSC', 'TFBS_ENCODE', 'TSS_Hoffman', 'H3K27ac_PGC2', 'H3K9ac_Trynka', 'Promoter_UCSC', 'H3K4me1_Trynka', 'H3K4me3_Trynka', 'Enhancer_Hoffman', 'Transcribed_Hoffman']
pos = hl.import_table("gs://qingbowang/ems_v1_test/{0}_positive_training_vg_beforeannot.tsv".format(tissue_name), impute=True).to_pandas()
neg = hl.import_table("gs://qingbowang/ems_v1_test/{0}_negative_training_vg_beforeannot.tsv".format(tissue_name), impute=True).to_pandas()

#keep the ids as v-g for final annotation
pos.index = pos.variant_id.str.replace("_b38","")
neg.index = neg.variant_id.str.replace("_b38","")
pos_vg_ix = pos.variant_id.str.replace("_b38","") + "_" + pos.gene_id
neg_vg_ix = neg.variant_id.str.replace("_b38","") + "_" + neg.gene_id

#and take the columns to use only
pos = pos[["tss_distance"] + cols_to_use_bin]
neg = neg[["tss_distance"] + cols_to_use_bin]
#and *1 to make it int instead of bool
pos = pos*1
neg = neg*1
#integrate basenji scores
bas_pos = []
bas_neg = []
for i in range(10+1):
    bas_pos_sub = hl.import_table("gs://qingbowang/ems_v1_test/{0}_sad_postraining_colchunk{1}.tsv".format(tissue_name, i), impute=True).to_pandas()
    bas_pos_sub.index = bas_pos_sub.ID
    del bas_pos_sub["ID"]
    bas_pos.append(bas_pos_sub)
    bas_neg_sub = hl.import_table("gs://qingbowang/ems_v1_test/{0}_sad_negtraining_colchunk{1}.tsv".format(tissue_name, i), impute=True).to_pandas()
    bas_neg_sub.index = bas_neg_sub.ID
    del bas_neg_sub["ID"]
    bas_neg.append(bas_neg_sub)
bas_pos = pd.concat(bas_pos, axis=1)
bas_neg = pd.concat(bas_neg, axis=1)

#bas has non-unique elements due to the way we subset the large table, but by definition, given the ID it is unique
#so we can remove duplicates:
bas_pos = bas_pos[~bas_pos.index.duplicated(keep='first')]
bas_neg = bas_neg[~bas_neg.index.duplicated(keep='first')]

#join them to pos and neg
pos = pos.join(bas_pos, how="left")
pos.index = pos_vg_ix #make the index v-g again (just to make it easier to track).
pos.dropna(subset=[bas_pos.columns[0]], inplace=True)#since sometimes NA due to throwing away hg19 no-match variants in gtex
pos.insert(0, "vg", pos.index) #and put the index as a column (otherwise hail drops it..)
neg = neg.join(bas_neg, how="left")
neg.index = neg_vg_ix #make the index v-g again (just to make it easier to track).
neg.dropna(subset=[bas_neg.columns[0]], inplace=True)#since sometimes NA due to throwing away hg19 no-match variants in gtex
neg.insert(0, "vg", neg.index)

hl.Table.from_pandas(pos).export("gs://qingbowang/ems_v1_test/{0}_positive_training_vg_annotated.tsv".format(tissue_name))
hl.Table.from_pandas(neg).export("gs://qingbowang/ems_v1_test/{0}_negative_training_vg_annotated.tsv".format(tissue_name))


#also annotate the binary features, tissue specific (and export as a separate file, because it is heavy)
pos = hl.import_table("gs://qingbowang/ems_v1_test/{0}_positive_training_vg_annotated.tsv".format(tissue_name)).select("vg")
neg = hl.import_table("gs://qingbowang/ems_v1_test/{0}_negative_training_vg_annotated.tsv".format(tissue_name)).select("vg")

pos = pos.add_index()
neg = neg.add_index() #just to match it later

pos = pos.annotate(ID = pos.vg.split("_")[0] + "_" + pos.vg.split("_")[1] + "_" + pos.vg.split("_")[2] + "_" + pos.vg.split("_")[3] + "_b38").key_by("ID") #stupid way but works..
neg = neg.annotate(ID = neg.vg.split("_")[0] + "_" + neg.vg.split("_")[1] + "_" + neg.vg.split("_")[2] + "_" + neg.vg.split("_")[3] + "_b38").key_by("ID")
marks = ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3" ,"H3K9ac", "H3K9me3"]
for m in marks:
    htsub = hl.read_table("gs://qingbowang/gtex_v8_wgs_binannot_pertissue_{0}.ht".format(m)).drop("locus")
    sub_pos = pos.join(htsub, how="left")
    sub_neg = neg.join(htsub, how="left")
    sub_pos = sub_pos.order_by("idx").drop("idx","ID") #match the original order
    sub_neg = sub_neg.order_by("idx").drop("idx", "ID")
    sub_pos.export("gs://qingbowang/ems_v1_test/{0}_positive_training_vg_roadmapannot_{1}.tsv".format(tissue_name, m))
    sub_neg.export("gs://qingbowang/ems_v1_test/{0}_negative_training_vg_roadmapannot_{1}.tsv".format(tissue_name, m))




