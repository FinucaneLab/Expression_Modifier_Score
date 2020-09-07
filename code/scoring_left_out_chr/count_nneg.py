# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')
import time as tm
import pandas as pd
tissue_name = "Whole_Blood"
chr_list = []
for i in range(1,22+1):
    chr_list.append("chr" + str(i))
chr_list.append("chrX")


nalls = []
nposs = []
nnegs = []
nnonposs = []
vg0 = hl.read_table("gs://qingbowang/ems_v1_test/{0}_allpairs_binannot_fmannot.ht".format(tissue_name))
vg0 = vg0.annotate(pos = ( (vg0.pp_susie>0.9) & (vg0.pp_fm>0.9)) )
vg0 = vg0.annotate(neg = ((vg0.pp_susie < 0.0001) & (vg0.pp_fm < 0.0001)) )
vg0 = vg0.annotate(chr=vg0.variant_id.split("_")[0])
for chr in chr_list:
    vg = vg0.filter(vg0.chr!=chr) #leaving out the chr
    N_all = vg.count()
    pos = vg.filter(vg.pos)
    n_pos = pos.count()
    neg = vg.filter((vg.neg))
    n_neg = neg.count()
    n_nonpos = N_all - n_pos
    nalls.append(N_all)
    nposs.append(n_pos)
    nnegs.append(n_neg)
    nnonposs.append(n_nonpos)
    print ("done {0}, {1}".format(chr, tm.ctime()))
df = pd.DataFrame({"n_all":nalls,"n_pos":nposs, "n_nonpos":nnonposs, "n_neg":nnegs})
df.index = chr_list
with hl.hadoop_open('gs://qingbowang/ems_v1_test/n_pos_and_neg_per_leftout_chr.tsv', 'w') as f: #index = chr to leave out
    df.to_csv(f, sep="\t")

