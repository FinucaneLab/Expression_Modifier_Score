# -*- coding: utf-8 -*-
__author__ = 'QingboWang'
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')

import pandas as pd
import numpy as np
import math

tissue_name = "Whole_Blood"



alphas = hl.import_table("gs://gtex_finemapping/v8/munged/GTEx_{0}_alphas.tsv.gz".format(tissue_name),  force_bgz=True, impute=True).repartition(400)
ems = hl.read_table("gs://qingbowang/ems_v1_test/ems_pcausal_gtexvg_all{0}.ht".format(tissue_name))

#for imputing, retrieve the p-random (first 10000 should be basically enough. theoretically this is constant other than precision error)
hthead = ems.head(10000)
hthead = hthead.annotate(p_causal_random = hthead.p_causal / hthead.confidence_gain)
prd = hthead.select("p_causal_random").to_pandas()
p_random = prd.p_causal_random.mean()#for imputing

#annotate ems for updating alpha
ems = ems.key_by("vg")
alphas = alphas.annotate(vg = alphas.variant_hg38 + "_" + alphas.gene).key_by("vg")
alphas = alphas.annotate(ems = ems[alphas.key].p_causal)
alphas = alphas.annotate(ems_was_missing = ~hl.is_defined(alphas.ems))
alphas = alphas.transmute(ems = hl.cond(alphas.ems_was_missing, p_random, alphas.ems))

#cap the ems so that min/max ratio = 100
M = alphas.aggregate(hl.agg.max(alphas.ems))
m = alphas.aggregate(hl.agg.min(alphas.ems))
import math
c = math.log(100, M/m)
alphas = alphas.annotate(ems = alphas.ems**c)

# the sum of alpha for each gene, which does not depend on the ems
sums_alpha = alphas.group_by("gene").aggregate(sum_alpha_1=hl.agg.sum(alphas.alpha1),
                                               sum_alpha_2=hl.agg.sum(alphas.alpha2),
                                               sum_alpha_3=hl.agg.sum(alphas.alpha3),
                                               sum_alpha_4=hl.agg.sum(alphas.alpha4),
                                               sum_alpha_5=hl.agg.sum(alphas.alpha5),
                                               sum_alpha_6=hl.agg.sum(alphas.alpha6),
                                               sum_alpha_7=hl.agg.sum(alphas.alpha7),
                                               sum_alpha_8=hl.agg.sum(alphas.alpha8),
                                               sum_alpha_9=hl.agg.sum(alphas.alpha9),
                                               sum_alpha_10=hl.agg.sum(alphas.alpha10))
alphas_orig = alphas #keeping the original version
alphas = alphas_orig #copy (deep)
colnames = list(pd.Series(np.arange(1,10+1)).apply(lambda x: "alpha"+str(x))) #so inefficient, but fine..
for i in range(len(colnames)):
    newname = "alpha_{0}_times_ems".format(i+1)
    alphas = alphas.annotate(tmp = alphas[colnames[i]] * alphas.ems)
    alphas = alphas.rename({"tmp": newname})  # and rename as new alpha

#prepare the sum for each gene (= each phenotype) to normalize by dividing by the sum
sums = alphas.group_by("gene").aggregate(sum_1=hl.agg.sum(alphas.alpha_1_times_ems),
                                                 sum_2=hl.agg.sum(alphas.alpha_2_times_ems),
                                                 sum_3=hl.agg.sum(alphas.alpha_3_times_ems),
                                                 sum_4=hl.agg.sum(alphas.alpha_4_times_ems),
                                                 sum_5=hl.agg.sum(alphas.alpha_5_times_ems),
                                                 sum_6=hl.agg.sum(alphas.alpha_6_times_ems),
                                                 sum_7=hl.agg.sum(alphas.alpha_7_times_ems),
                                                 sum_8=hl.agg.sum(alphas.alpha_8_times_ems),
                                                 sum_9=hl.agg.sum(alphas.alpha_9_times_ems),
                                                 sum_10=hl.agg.sum(alphas.alpha_10_times_ems) #doesn't scale .. but fine for just 10 effects
                                                 )

alphas = alphas.key_by("gene").join(sums, how="left")
alphas = alphas.join(sums_alpha, how="left")
#and divide by the sum, and multiply by the sum alpha, since we don't want to screw up when the sum of alpha is not 1
for i in range(len(colnames)):
    newname ="alpha_updated_{0}".format(i+1)
    alphas = alphas.annotate(tmp = alphas["alpha_{0}_times_ems".format(i+1)] / alphas["sum_{0}".format(i+1)] * alphas["sum_alpha_{0}".format(i+1)])
    alphas = alphas.rename({"tmp": newname})  # and rename as new alpha
#we just want the sum of unif_alpha to be the sum of updated_alpha.

#then get the final updated pip by 1 - \prod{alphas}
#-> where, alpha = alpha_u when no signal, alpha_f only if the alpha (the csid) is meaningful



#to do so, first get the cs list
csid = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/Whole_Blood_csid_list_per_gene.ht")
alphas = alphas.key_by("gene").join(csid, how="left")

for i in range(len(colnames)):
    newname = "alpha_updated_{0}".format(i + 1)
    alphas = alphas.annotate(tmp = hl.cond( alphas["csid_list"].contains(hl.int64(i+1)), alphas[newname], alphas["alpha{0}".format(i+1)]))
    alphas = alphas.drop(newname)
    alphas = alphas.rename({"tmp": newname})  # and rename as new alpha (this could be alpha_u or alpha_f depending on the cs, but mostly alpha_u this is)


alphas = alphas.annotate(prod = 1)
for i in range(len(colnames)):
    newname = "alpha_updated_{0}".format(i+1)
    alphas = alphas.annotate(prod = alphas.prod * (1-alphas[newname])) #multiplying one by one
alphas = alphas.annotate(updated_pip = 1 - alphas.prod)
#and finally write as a hail table
alphas.write("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_alpha_filtered_by_cs_minmax100.ht".format(tissue_name), overwrite=True)

#sum stats
ht = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_alpha_filtered_by_cs_minmax100.ht".format(tissue_name))

ht = ht.annotate(pip_bin_unif_prior = hl.case()
                                             .when(ht.pip > 0.9, 9)
                                             .when(ht.pip > 0.5, 5)
                                             .when(ht.pip > 0.1, 1)
                                             .when(ht.pip > 0.01, 0)
                                             .default(-1),
                 pip_bin_ems_as_a_prior = hl.case()
                 .when(ht.updated_pip > 0.9, 9)
                 .when(ht.updated_pip > 0.5, 5)
                 .when(ht.updated_pip > 0.1, 1)
                 .when(ht.updated_pip > 0.01, 0)
                 .default(-1)
                 )
ht.group_by("pip_bin_unif_prior", "pip_bin_ems_as_a_prior").aggregate(n=hl.agg.count()).export("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_overview_cutoff_cs_minmax100.tsv".format(tissue_name))


#also vs pip_f
ht = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_alpha_filtered_by_cs_minmax100.ht".format(tissue_name))

ht1 = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_alpha_filtered_by_cs.ht".format(tissue_name)).key_by("vg")

ht = ht.annotate(pip_bin_fullems_prior = hl.case()
                                             .when(ht1[ht.vg].updated_pip > 0.9, 9)
                                             .when(ht1[ht.vg].updated_pip > 0.5, 5)
                                             .when(ht1[ht.vg].updated_pip > 0.1, 1)
                                             .when(ht1[ht.vg].updated_pip > 0.01, 0)
                                             .default(-1),
                 pip_bin_ems_as_a_prior = hl.case()
                 .when(ht.updated_pip > 0.9, 9)
                 .when(ht.updated_pip > 0.5, 5)
                 .when(ht.updated_pip > 0.1, 1)
                 .when(ht.updated_pip > 0.01, 0)
                 .default(-1)
                 )
ht.group_by("pip_bin_fullems_prior", "pip_bin_ems_as_a_prior").aggregate(n=hl.agg.count()).export("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_overview_cutoff_cs_minmax100_vsfull.tsv".format(tissue_name))
