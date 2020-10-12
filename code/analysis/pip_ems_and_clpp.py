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


#output the ones with non-trivial signal for PIP_unif or PIP_EMS
for tissue_name in tissues:
    ht = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_alpha_filtered_by_cs.ht".format(tissue_name))
    ht.filter((ht.pip>0.1) | (ht.updated_pip>0.1)).write("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_thres01.ht".format(tissue_name), overwrite=True)

## PIP_EMS vs PIP_unif stats

#get the list of newly identified putative causal eQTLs
dfall = []
for tissue_name in tissues:
    df = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_thres01.ht".format(tissue_name))
    df = df.filter((df.pip<0.9) & (df.updated_pip>0.9)).to_pandas()
    df["tissue"] = tissue_name
    dfall.append(df)
dfall = pd.concat(dfall, axis=0)
fn = "gs://qingbowang/tmp/newly_found_putative_causals.tsv"
with hl.hadoop_open(fn, 'w') as f:
    dfall.to_csv(f, sep="\t")

#also lost putative causal eQTLs:
dfall = []
for tissue_name in tissues:
    df = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_thres01.ht".format(tissue_name))
    df = df.filter((df.pip>0.9) & (df.updated_pip<0.9)).to_pandas()
    df["tissue"] = tissue_name
    dfall.append(df)
dfall = pd.concat(dfall, axis=0)
fn = "gs://qingbowang/tmp/lost_putative_causals.tsv"
with hl.hadoop_open(fn, 'w') as f:
    dfall.to_csv(f, sep="\t")
#also consistent ones
dfall = []
for tissue_name in tissues:
    df = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_thres01.ht".format(tissue_name))
    df = df.filter((df.pip>0.9) & (df.updated_pip>0.9)).to_pandas()
    df["tissue"] = tissue_name
    dfall.append(df)
dfall = pd.concat(dfall, axis=0)
fn = "gs://qingbowang/tmp/consistent_putative_causals.tsv"
with hl.hadoop_open(fn, 'w') as f:
    dfall.to_csv(f, sep="\t")

#stats:
fn = "gs://qingbowang/tmp/newly_found_putative_causals.tsv"
with hl.hadoop_open(fn, 'r') as f:
    df = pd.read_csv(f, sep="\t")
fn = "gs://qingbowang/tmp/newly_found_putative_causals_n.tsv"
with hl.hadoop_open(fn, 'w') as f:
    df.tissue.value_counts().sort_values().to_csv(f, sep="\t")


##co-localization analysis stats:

#get the comp trait signals:
comp = hl.read_table("gs://qingbowang/UKBB_94traits_pp.ht")
comp.head(5).show()
comp.filter((comp.pip>0.1)&(comp.method=="SUSIE")).write("gs://qingbowang/UKBB_pp_susie_01.ht", overwrite=True)


#output the clpp, for non-trivial ones:
comp = hl.read_table("gs://qingbowang/UKBB_pp_susie_01.ht")
comp = comp.key_by().select("trait","pip", "hg38_ID").rename({"pip":"pip_comp"}).to_pandas()
comp.index = comp.hg38_ID
dfall = []
for tissue_name in tissues:
    eq = hl.read_table("gs://qingbowang/ems_v1_test/updated_pips/{0}_updated_pips_thres01.ht".format(tissue_name))
    eq = eq.key_by().select("gene","variant_hg38","pip","updated_pip").rename({"pip":"pip_u", "updated_pip":"pip_ems"}).to_pandas()
    eq.index = eq.variant_hg38.str.replace("_b38","")
    eq["tissue"] = tissue_name
    intr =comp.join(eq, how="inner")
    if tissue_name == tissues[0]:
        print(intr.head(10))
    dfall.append(intr)
dfall = pd.concat(dfall, axis=0)
dfall["clpp"] = dfall.pip_comp*dfall.pip_u
dfall["ems_clpp"] = dfall.pip_comp*dfall.pip_ems
fn = "gs://qingbowang/tmp/table_for_clpp.tsv"
with hl.hadoop_open(fn, 'w') as f:
    dfall[(dfall.clpp>0.1)|(dfall.ems_clpp>0.1)].to_csv(f, sep="\t")


#final num. pipu and pipems
#(in local):
df = pd.read_csv("~/Downloads/table_for_clpp.tsv", sep="\t")
print (df[df.clpp>0.1].gene.value_counts().shape)
print (df[df.ems_clpp>0.1].gene.value_counts().shape)

print ((df.clpp>0.1).sum())
print ((df.ems_clpp>0.1).sum())
print ((df.ems_clpp>0.1).sum() - (df.clpp>0.1).sum())

#gene-trait:
print (df[df.clpp>0.1].groupby(["gene","trait"]).size().shape)
print (df[df.ems_clpp>0.1].groupby(["gene","trait"]).size().shape)

#also dapg CLPP
dp = hl.import_table("gs://qingbowang/GTEx_v8_finemapping_DAPG.txt.gz", impute=True, force=True)
dp.repartition(400).write("gs://qingbowang/tmp/dapg_alltissue.ht")

dp = hl.read_table("gs://qingbowang/tmp/dapg_alltissue.ht")
dp.filter(dp.variant_pip>0.1).write("gs://qingbowang/tmp/dapg_alltissue_pip01.ht")

dp = dp.select("variant_id","variant_pip", "tissue_id", "gene_id").rename({"variant_pip":"pip_dapg"}).to_pandas()
dp.index = dp.variant_id.str.replace("_b38","")
comp = hl.read_table("gs://qingbowang/UKBB_pp_susie_01.ht")
comp = comp.key_by().select("trait","pip", "hg38_ID").rename({"pip":"pip_comp"}).to_pandas()
comp.index = comp.hg38_ID
dfall =comp.join(dp, how="inner")
dfall["dapg_clpp"] = dfall.pip_comp*dfall.pip_dapg
fn = "gs://qingbowang/tmp/table_for_clpp_dapg.tsv"
with hl.hadoop_open(fn, 'w') as f:
    dfall[(dfall.dapg_clpp>0.1)].to_csv(f, sep="\t")



#final stats (in local):
df = pd.read_csv("~/Downloads/table_for_clpp_dapg.tsv", sep="\t")
print (df[df.dapg_clpp>0.1].gene_id.value_counts().shape)
print ((df.dapg_clpp>0.1).sum())
print (df[df.dapg_clpp>0.1].groupby(["gene_id","trait"]).size().shape)
