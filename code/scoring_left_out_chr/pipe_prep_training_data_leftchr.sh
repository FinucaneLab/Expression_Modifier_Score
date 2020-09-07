#!/bin/bash


#specify the tissue
tissue_name=Whole_Blood
leftout_chr=$1
chr_name_for_cluster=$(tr '[A-Z]' '[a-z]' <<< $leftout_chr)

#0.1. spin up a cluster
hailctl dataproc start qbwcluster-$chr_name_for_cluster -w 2 --max-idle 10m --pkgs matplotlib,seaborn,sklearn ; #small is fine

#3. do the training

#3.1. select subset of tissue-specific binary features by marginal f1 value -> in hail!
hailctl dataproc submit qbwcluster-$chr_name_for_cluster select_tissue_spec_bin_feat.py Whole_Blood $leftout_chr;

#3.2. select subset of basenji features by marginal ROC
hailctl dataproc submit qbwcluster-$chr_name_for_cluster select_bas_feat.py Whole_Blood $leftout_chr ;

