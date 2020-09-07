# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

#downstream = turn the score to "p-causal"

import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')
import pandas as pd
import time as tm

tissue_name = "Whole_Blood"

chr_list = []
for i in range(1,22+1):
    chr_list.append("chr" + str(i))
chr_list.append("chrX")

for chr in chr_list:
    print ("starting {0}, {1}".format(chr, tm.ctime()))

    fn = "gs://qingbowang/ems_v1_test/ems_p_causal_interpolated_{0}_leftout_{1}.tsv".format(tissue_name, chr)
    with hl.hadoop_open(fn, 'r') as f:
        pcausal = pd.read_csv(f, sep="\t", index_col=0)
    pcausal["rf_score_bin"] = pcausal.index
    del pcausal["rf_score_bin.1"] #duplicated columns
    
    pcausal = hl.Table.from_pandas(pcausal)
    pcausal = pcausal.transmute(rf_score_bin = hl.format('%.3f', pcausal["rf_score_bin"]))
    
    #score the chunk at once:
    df = hl.import_table("gs://qingbowang/ems_v1_test/ems_rawscore_gtexvg_all_Whole_Blood_{0}.tsv.gz".format(chr), force=True, impute=True)
    df = df.repartition(80)
    df = df.annotate(rf_score_bin = hl.format('%.3f', df["0"]))
    pcausal = pcausal.key_by("rf_score_bin")
    df = df.key_by("rf_score_bin")
    df = df.join(pcausal, how="left")
    df = df.rename({'' : 'vg', '0' : 'rf_score_raw'})
    df.write("gs://qingbowang/ems_v1_test/tmp/ems_pcausal_gtexvg_all_Whole_Blood_chr{0}.ht".format(chr), overwrite=True)


#concat all the chrom
dfall = []
for chr in chr_list:
    df = hl.read_table("gs://qingbowang/ems_v1_test/tmp/ems_pcausal_gtexvg_all_Whole_Blood_chr{0}.ht".format(chr))
    dfall.append(df)
#and finally, combine them all, and write
dfall = dfall[0].union(*dfall[1:])
dfall.write("gs://qingbowang/ems_v1_test/ems_pcausal_gtexvg_all_Whole_Blood_left_each_chr.ht", overwrite=True)






