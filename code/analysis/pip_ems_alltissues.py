# -*- coding: utf-8 -*-
__author__ = 'QingboWang'
import hail as hl
import hail.expr.aggregators as agg
hl.init(default_reference='GRCh38')

import pandas as pd
import numpy as np

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
            "Testis",
            "Artery_Aorta",
            "Brain_Amygdala",
            "Brain_Anterior_cingulate_cortex_BA24",
            "Brain_Caudate_basal_ganglia",
            "Brain_Cerebellar_Hemisphere",
            "Brain_Frontal_Cortex_BA9",
            "Brain_Hippocampus",
            "Brain_Hypothalamus",
            "Brain_Putamen_basal_ganglia",
            "Brain_Spinal_cord_cervical_c-1",
            "Brain_Substantia_nigra",
            "Cells_Cultured_fibroblasts",
            "Cells_EBV-transformed_lymphocytes",
            "Kidney_Cortex",
            "Minor_Salivary_Gland",
            "Ovary",
            "Small_Intestine_Terminal_Ileum",
            "Uterus",
            "Vagina"]#all 49 tissues
for tissue_name in tissues:

    #get the alpha for each tissue
    alphas = hl.import_table("gs://gtex_finemapping/v8/munged/GTEx_{0}_alphas.tsv.gz".format(tissue_name),  force_bgz=True, impute=True).repartition(400)
    alphas.group_by("gene").aggregate(n = hl.agg.count()).export("gs://qingbowang/ems_v1_test/updated_pips/{0}_n_vars_per_gene_alpha.tsv".format(tissue_name))

    fn = "gs://qingbowang/ems_v1_test/updated_pips/{0}_n_vars_per_gene_alpha.tsv".format(tissue_name)
    with hl.hadoop_open(fn, 'r') as f:
        ng = pd.read_csv(f, sep="\t")
    ng = hl.Table.from_pandas(ng).key_by("gene")
    alphas = hl.import_table("gs://gtex_finemapping/v8/munged/GTEx_{0}_alphas.tsv.gz".format(tissue_name),  force_bgz=True, impute=True).repartition(400)
    alphas.key_by("gene").join(ng, how="left").write("gs://qingbowang/ems_v1_test/updated_pips/GTEx_{0}_alphas_and_n_vars.ht".format(tissue_name), overwrite=True)
    ht = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/GTEx_{0}_alphas_and_n_vars.ht".format(tissue_name))

    #for each gene, get the meaningful csid:
    ht.group_by("gene","cs_id").aggregate(n = hl.agg.count()).export("gs://qingbowang/ems_v1_test/updated_pips/{0}_csid_cnt_per_gene.tsv".format(tissue_name))
    fn = "gs://qingbowang/ems_v1_test/updated_pips/{0}_csid_cnt_per_gene.tsv".format(tissue_name)
    with hl.hadoop_open(fn, 'r') as f:
        df = pd.read_csv(f, sep="\t")
    l = df.groupby('gene')['cs_id'].apply(list)
    l = pd.DataFrame(l)
    l["gene"] = l.index
    l = l.rename(columns={'cs_id': 'csid_list'})
    ht = hl.Table.from_pandas(l).key_by("gene")
    ht.write("gs://qingbowang/ems_v1_test/updated_pips/{0}_csid_list_per_gene.ht".format(tissue_name), overwrite=True) #row = gene, column = csid list
  

    #and use that for pip_f update:
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
                                                     sum_10=hl.agg.sum(alphas.alpha_10_times_ems)
                                                     )

    alphas = alphas.key_by("gene").join(sums, how="left")
    alphas = alphas.join(sums_alpha, how="left")
    #and divide by the sum, and multiply by the sum alpha
    for i in range(len(colnames)):
        newname ="alpha_updated_{0}".format(i+1)
        alphas = alphas.annotate(tmp = alphas["alpha_{0}_times_ems".format(i+1)] / alphas["sum_{0}".format(i+1)] * alphas["sum_alpha_{0}".format(i+1)])
        alphas = alphas.rename({"tmp": newname})  # and rename as new alpha

    #then get the final updated pip by 1 - \prod{alphas}
    #where, alpha = alpha_u when no signal, alpha_f only if the alpha (the csid) is meaningful

    #to do so, first get the cs list
    csid = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/{0}_csid_list_per_gene.ht".format(tissue_name))
    alphas = alphas.key_by("gene").join(csid, how="left")

    for i in range(len(colnames)):
        newname = "alpha_updated_{0}".format(i + 1)
        alphas = alphas.annotate(tmp = hl.cond( alphas["csid_list"].contains(hl.int64(i+1)), alphas[newname], alphas["alpha{0}".format(i+1)]))
        alphas = alphas.drop(newname)
        alphas = alphas.rename({"tmp": newname})  # and rename as new alpha (this could be alpha_u or alpha_f depending on the cs, but mostly alpha_u this is)
    #do the \prod
    alphas = alphas.annotate(prod = 1)
    for i in range(len(colnames)):
        newname = "alpha_updated_{0}".format(i+1)
        alphas = alphas.annotate(prod = alphas.prod * (1-alphas[newname])) #multiplying one by one
    alphas = alphas.annotate(updated_pip = 1 - alphas.prod)
    #and finally write as a hail table
    alphas.select("variant_hg38","variant","vg","pip", "ems", "updated_pip").write("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_alpha_filtered_by_cs.ht".format(tissue_name), overwrite=True)

