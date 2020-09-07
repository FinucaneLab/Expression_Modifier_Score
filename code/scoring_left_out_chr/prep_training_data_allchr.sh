#!/bin/bash

(bash ./pipe_prep_training_data_leftchr.sh chr1) &
(bash ./pipe_prep_training_data_leftchr.sh chr2) &
(bash ./pipe_prep_training_data_leftchr.sh chr3) &
(bash ./pipe_prep_training_data_leftchr.sh chr4) &
(bash ./pipe_prep_training_data_leftchr.sh chr5) &
(bash ./pipe_prep_training_data_leftchr.sh chr6) &
(bash ./pipe_prep_training_data_leftchr.sh chr7) &
(bash ./pipe_prep_training_data_leftchr.sh chr8) &
(bash ./pipe_prep_training_data_leftchr.sh chr9) &
(bash ./pipe_prep_training_data_leftchr.sh chr10) &
(bash ./pipe_prep_training_data_leftchr.sh chr11) &
(bash ./pipe_prep_training_data_leftchr.sh chr12) &
(bash ./pipe_prep_training_data_leftchr.sh chr13) &
(bash ./pipe_prep_training_data_leftchr.sh chr14) &
(bash ./pipe_prep_training_data_leftchr.sh chr15) &
(bash ./pipe_prep_training_data_leftchr.sh chr16) &
(bash ./pipe_prep_training_data_leftchr.sh chr17) &
(bash ./pipe_prep_training_data_leftchr.sh chr18) &
(bash ./pipe_prep_training_data_leftchr.sh chr19) &
(bash ./pipe_prep_training_data_leftchr.sh chr20) &
(bash ./pipe_prep_training_data_leftchr.sh chr21) &
(bash ./pipe_prep_training_data_leftchr.sh chr22) &
(bash ./pipe_prep_training_data_leftchr.sh chrX) &
wait


#do the downstream in a small cluster
hailctl dataproc start qbwcluster-matome -w 2 --max-idle 10m --pkgs matplotlib
hailctl dataproc submit qbwcluster-matome ./select_tissue_spec_bin_feat_matome.py
hailctl dataproc submit qbwcluster-matome ./select_bas_feat_matome.py
hailctl dataproc stop qbwcluster-matome

#count the nneg here:
hailctl dataproc start qbwcluster -w 20 --max-idle 10m
hailctl dataproc submit qbwcluster ./count_nneg.py
hailctl dataproc stop qbwcluster
