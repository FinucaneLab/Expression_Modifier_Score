# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')

#non-tissue specific. already done
#0.1. add binary annotations
def annotate_sites(ht, bed_dir = "gs://gnomad-qingbowang/finucane_et_al_hg38_ht/"): #currently only supports hg38. Returns annotated hail table
    #per functional annotation, annotate binary
    ls = hl.hadoop_ls(bed_dir)
    for i in range(len(ls)):
        func_interval = hl.read_table(ls[i]["path"])
        func_name = ls[i]["path"].split("/")[-1].split(".")[0]
        ht = ht.annotate(func_name_tmp = hl.is_defined(func_interval[ht.locus])) #temp for the column name
        ht = ht.rename({"func_name_tmp":func_name}) #and change it to the functional anntation name itself
        print ("done {0}".format(func_name))
    return (ht)
v = hl.read_table("gs://qingbowang/gtex_v8_wgs.ht")
v = v.key_by(v.variant_id).select()
v = v.annotate(locus=hl.locus( v.variant_id.split("_")[0], hl.int32(v.variant_id.split("_")[1])))
vsub = annotate_sites(v)
vsub.write("gs://qingbowang/gtex_v8_wgs_binannot.ht", overwrite=True)


#tissue specific annotations as well
def annotate_sites_specific_histonemark(ht, bed_dir = "gs://gnomad-qingbowang/finucane_et_al_hg38_ht/", mark="H3K4me1"):
    #per histone mark, since otherwise this will blow up
    #per functional annotation, annotate binary
    ls = hl.hadoop_ls(bed_dir)
    for i in range(len(ls)):
        func_interval = hl.read_table(ls[i]["path"])
        func_name = ls[i]["path"].split("/")[-1].replace("_narrowpeak.ht", "")
        his_name = func_name.split("-")[-1]
        if his_name ==mark:
            ht = ht.annotate(func_name_tmp = hl.is_defined(func_interval[ht.locus])) #temp for the column name
            ht = ht.rename({"func_name_tmp":func_name}) #and change it to the functional anntation name itself
            print ("done {0}".format(func_name))
    return (ht)
# ~= 127*7 ~= 900 columns at once? Yeah let's do.
#no, one histone mark at a time
#for h in ["H3K9me3", "H3K36me3", "H3K27me3", "H3K4me3", "H3K4me1", "H3K27ac", "H3K9ac"]: #order confusing, make it uniform across pipeline
for h in ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3"]:
    vsub = annotate_sites_specific_histonemark(v, bed_dir="gs://qingbowang/roadmap_bed_narrowpeak_hg38_ht", mark=h)
    vsub.write("gs://qingbowang/gtex_v8_wgs_binannot_pertissue_{0}.ht".format(h), overwrite=True)
