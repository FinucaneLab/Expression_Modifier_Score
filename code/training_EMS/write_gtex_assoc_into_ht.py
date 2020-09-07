# -*- coding: utf-8 -*-
__author__ = 'QingboWang'


import sys
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')
tissue_name = sys.argv[1]
t = hl.import_table("gs://gtex_finemapping/v8/GTEx_Analysis_v8_eQTL_all_associations/{0}.allpairs.txt.gz".format(tissue_name),
                    force=True, types={'tss_distance': hl.tint32, 'maf': hl.tfloat64, "pval_nominal":hl.tfloat64, "slope":hl.tfloat64})
t = t.repartition(800).write("gs://qingbowang/ems_v1_test/{0}_allpairs.ht".format(tissue_name), overwrite=True)

