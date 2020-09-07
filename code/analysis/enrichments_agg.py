import sys
import pandas as pd
import numpy as np
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')



#enrichment of baseline features:
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
            "Testis"]
his_to_use = ["H3K27ac_PGC2","H3K9ac_Trynka","H3K4me1_Trynka","H3K4me3_Trynka"]
nc_to_use = ["DHS_Trynka","TFBS_ENCODE","Promoter_UCSC","Enhancer_Hoffman"]
c_to_use = ["UTR_5_UCSC","UTR_3_UCSC", "Coding_UCSC","Intron_UCSC"]

for tissue_name in tissues:
    vg = hl.read_table("gs://qingbowang/ems_v1_test/{0}_allpairs_binannot_fmannot.ht".format(tissue_name))
    vg = vg.annotate(pip_bin = hl.case()
                                 .when(vg.pp_susie>0.9, 9)
                                 .when(vg.pp_susie>0.5, 5)
                                 .when(vg.pp_susie>0.1, 1)
                                 .when(vg.pp_susie>0.01, 0)
                                 .default(-1) )
    for h in his_to_use + nc_to_use + c_to_use:
        df = vg.group_by(h, "pip_bin").aggregate(n=hl.agg.count())
        df.export("gs://qingbowang/ems_v1_test/agg_for_fig1/{0}_{1}_agg.tsv".format(tissue_name, h))


#downstream: collecting them as a single dataframe
import pandas as pd
import numpy as np
for tissue_name in tissues:
    for h in his_to_use + nc_to_use + c_to_use:
        with hl.hadoop_open("gs://qingbowang/ems_v1_test/agg_for_fig1/{0}_{1}_agg.tsv".format(tissue_name, h), 'r') as f:
            t = pd.read_csv(f, sep="\t")
        if h == his_to_use[0]:#if the first, need to annotate the "total"
            df = pd.DataFrame(t.groupby("pip_bin").n.sum())  # total, before annotating things
            df.columns = ["n_total"]
            df["frac_naive"] = df.n_total / df.n_total.sum()
            dfn = df.copy(deep=True) #to output the ns
            dferr = df.copy(deep=True) #to output the err
        dfpos = t[t[h]]  # positive for that feature
        dfpos.index = dfpos.pip_bin
        dfpos["n_pos"] = dfpos.n
        dfpos["frac_cond"] = dfpos.n_pos / dfpos.n_pos.sum() #fraction conditioning on that histone mark
        dfpos["err_cond"] = np.sqrt(dfpos.frac_cond * (1 - dfpos.frac_cond) / dfpos.n_pos.sum())
        dfpos["enr_cond"] = dfpos.frac_cond / df.frac_naive #enrichment
        dfpos["enr_cond_err"] = dfpos.err_cond / df.frac_naive  # although way of putting this error is questionable
        df["enr_{0}".format(h)] = dfpos["enr_cond"]
        dfn["n_{0}".format(h)] = dfpos["n_pos"]
        dferr["err_{0}".format(h)] = dfpos["enr_cond_err"]
    #and save
    with hl.hadoop_open("gs://qingbowang/ems_v1_test/agg_for_fig1/{0}_binfeatures_enrichments.tsv".format(tissue_name), 'w') as f:
        df.to_csv(f, sep="\t")
    with hl.hadoop_open("gs://qingbowang/ems_v1_test/agg_for_fig1/{0}_binfeatures_ns.tsv".format(tissue_name), 'w') as f:
        dfn.to_csv(f, sep="\t")
    with hl.hadoop_open("gs://qingbowang/ems_v1_test/agg_for_fig1/{0}_binfeatures_errorbars.tsv".format(tissue_name), 'w') as f:
        dferr.to_csv(f, sep="\t")


#TSS enrichment:
tissue_name = "Whole_Blood"
vg = hl.read_table("gs://qingbowang/ems_v1_test/{0}_allpairs_binannot_fmannot.ht".format(tissue_name))
vg = vg.annotate(abs_tss_dist = hl.abs(vg.tss_distance))
vg = vg.annotate(pip_bin = hl.case()
                                 .when(vg.pp_susie>0.9, 9)
                                 .when(vg.pp_susie>0.5, 5)
                                 .when(vg.pp_susie>0.1, 1)
                                 .when(vg.pp_susie>0.01, 0)
                                 .default(-1),
                 tss_bin = hl.case()
                                    .when(vg.abs_tss_dist<10**2,2)#100
                                    .when(vg.abs_tss_dist<10**2.5,2.5)
                                    .when(vg.abs_tss_dist<10**3,3) #1k
                                    .when(vg.abs_tss_dist<10**3.5,3.5)
                                    .when(vg.abs_tss_dist<10**4,4)#10k
                                    .when(vg.abs_tss_dist<10**4.5,4.5)
                                    .when(vg.abs_tss_dist<10**5,5)#100k
                                    .when(vg.abs_tss_dist<10**5.5,5.5)#100k
                                    .default(6) #~1m
                 )
df = vg.group_by("tss_bin", "pip_bin").aggregate(n=hl.agg.count())
df.export("gs://qingbowang/ems_v1_test/agg_for_fig1/{0}_tss_dist_agg.tsv".format(tissue_name))

#Basenji feature enrichment:
bins = np.logspace(0,2,9)
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
i = 0
for f in fs:
    htsub = hl.read_table(f)
    rowstouse = np.intersect1d(list(htsub.row), bas_for_plot)
    htsub = htsub.select(*rowstouse)
    htsub.write("gs://qingbowang/ems_v1_test/sad_for_plot_sub{0}.ht".format(i), overwrite=True)
    htsub = hl.read_table("gs://qingbowang/ems_v1_test/sad_for_plot_sub{0}.ht".format(i))
    if i==0:
        ht = htsub
    else:
        ht = ht.join(htsub, how="left")
    i += 1
import numpy as np
for tissue_name in ["Whole_Blood"]:
    vg = hl.read_table("gs://qingbowang/ems_v1_test/{0}_allpairs_binannot_fmannot.ht".format(tissue_name)).key_by("variant_id")
    N = vg.count()
    # and aggregate
    vg = vg.annotate(pip_bin=hl.case()
                     .when(vg.pp_susie > 0.9, 9)
                     .when(vg.pp_susie > 0.5, 5)
                     .when(vg.pp_susie > 0.1, 1)
                     .when(vg.pp_susie > 0.01, 0)
                     .default(-1))
    for j in range(i):
        ht = hl.read_table("gs://qingbowang/ems_v1_test/sad_for_plot_sub{0}.ht".format(j))
        basall = list(ht.row)[1:]  # not including the "variant ID"
        if "CNhs10862" in basall:#this is the one to plot
            h = "CNhs10862"
            ht = ht.annotate(ID_b38 = ht.ID + "_b38")
            ht = ht.key_by("ID_b38")
            vgj = vg.join(ht, how="left")
            vgj = vgj.annotate(abs_score = hl.abs(vgj[h]))
            vgj = vgj.annotate(bas_bin = hl.case()
                                               .when(vgj.abs_score > bins[-1], 9)
                                               .when(vgj.abs_score > bins[-2], 8)
                                               .when(vgj.abs_score > bins[-3], 7)
                                               .when(vgj.abs_score > bins[-4], 6)
                                               .when(vgj.abs_score > bins[-5], 5)
                                               .when(vgj.abs_score > bins[-6], 4)
                                               .when(vgj.abs_score > bins[-7], 3)
                                               .when(vgj.abs_score > bins[-8], 2)
                                               .when(vgj.abs_score > bins[-9], 1)
                                               .default(0)
                                   )
            df = vgj.group_by("bas_bin", "pip_bin").aggregate(n=hl.agg.count())
            df.export("gs://qingbowang/ems_v1_test/agg_for_fig1/{0}_{1}_basenji_agg_forplot.tsv".format(tissue_name, h))



#quantitative enrichment of pip itself (n: log)
tissue_name = "Whole_Blood"
ems = hl.read_table("gs://qingbowang/ems_v1_test/ems_pcausal_gtexvg_all{0}.ht".format(tissue_name)).key_by("vg")
vg = hl.read_table("gs://qingbowang/ems_v1_test/{0}_allpairs_binannot_fmannot.ht".format(tissue_name))
vg = vg.annotate(vg = vg.variant_id+"_"+vg.gene_id)
vg = vg.annotate(ems_bin = hl.ceil(ems[vg.vg].confidence_gain_log10))
vg = vg.annotate(pip_bin = hl.case()
                                 .when(vg.pp_susie>0.9, 9)
                                 .when(vg.pp_susie>0.5, 5)
                                 .when(vg.pp_susie>0.1, 1)
                                 .when(vg.pp_susie>0.01, 0)
                                 .default(-1)
                 )
df = vg.group_by("pip_bin", "ems_bin").aggregate(n=hl.agg.count())
df.export("gs://qingbowang/ems_v1_test/agg_for_fig1/{0}_ems_agg.tsv".format(tissue_name))
