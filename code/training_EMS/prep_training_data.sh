#!/bin/bash


#0.1. add binary annotations to the original variants file
python3 bin_annotate_gtex_var.py

#specify the tissue
tissue_name=$1
# _ and uppercase cannot be used for cluster name, so replace
tissue_name_for_cluster="${tissue_name//_/}"
tissue_name_for_cluster=$(tr '[A-Z]' '[a-z]' <<< $tissue_name_for_cluster)

#0.1. spin up a cluster
hailctl dataproc start qbwcluster-$tissue_name_for_cluster -w 2 --max-idle 10m ;

#0.2. write the files as ht
hailctl dataproc submit qbwcluster-$tissue_name_for_cluster write_gtex_assoc_into_ht.py $tissue_name ;

#1. annotate the baseline annotations and pip
hailctl dataproc submit qbwcluster-$tissue_name_for_cluster annotate_bin_and_pp.py $tissue_name ;

#2. define the training data, annotate basenji scores, and export
hailctl dataproc submit qbwcluster-$tissue_name_for_cluster export_pos_and_neg_for_training.py $tissue_name ;

#3. do the training

#3.1. select subset of tissue-specific binary features by marginal f1 value
hailctl dataproc submit qbwcluster-$tissue_name_for_cluster select_tissue_spec_bin_feat.py $tissue_name ;

#3.2. select subset of basenji features by marginal ROC
hailctl dataproc submit qbwcluster-$tissue_name_for_cluster select_bas_feat.py $tissue_name ;

#3.3. combine those features to use, and further perform final best-tuned algorithm
#also output the predictor itself / ems score conversion thing here as well

#4. do the scoring of the rest of GTEx tissues



