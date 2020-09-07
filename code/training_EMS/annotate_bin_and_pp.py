# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

import sys
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')

tissue_name = sys.argv[1]

#1.1 load all the GTEx variants (non tissue specific)
v = hl.read_table("gs://qingbowang/gtex_v8_wgs_binannot.ht")
v = v.key_by("variant_id")

#1.2. load all the gene-variant pairs for your tissue of interest
vg = hl.read_table("gs://qingbowang/ems_v1_test/{0}_allpairs.ht".format(tissue_name))

#2. annotate binary features, for each v-g pair
vg = vg.key_by("variant_id").join(v, how="left")
vg.write("gs://qingbowang/ems_v1_test/{0}_allpairs_binannot.ht".format(tissue_name), overwrite=True)

#2.2. also fine mapping prob
vg = hl.read_table("gs://qingbowang/ems_v1_test/{0}_allpairs_binannot.ht".format(tissue_name)).key_by("variant_id","gene_id")

#annotate fm prob

#read the fm file
fm = hl.read_table("gs://qingbowang/susie_results/GTEx_30tissues_pip_FINEMAP.ht").key_by("variant_hg38","gene")
sus = hl.read_table("gs://qingbowang/susie_results/GTEx_30tissues_pip_SuSiE.ht").key_by("variant_hg38","gene")

#filter to the tissue of interest
sus = sus.filter(sus.tissue==tissue_name)
fm = fm.filter(fm.tissue==tissue_name)

#annotate the v-g
vg = vg.annotate(pp_susie = hl.cond(hl.is_defined(sus[vg.key].pip), sus[vg.key].pip, 0),
                 pp_fm = hl.cond(hl.is_defined(fm[vg.key].pip), fm[vg.key].pip, 0))
#and export
vg.write("gs://qingbowang/ems_v1_test/{0}_allpairs_binannot_fmannot.ht".format(tissue_name), overwrite=True)







