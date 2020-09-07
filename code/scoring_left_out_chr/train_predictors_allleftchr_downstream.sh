#!/usr/bin/env bash

#gcloud beta compute らへんは train_predictors_alltissues.sh と同じ

#localで:
gsutil cp ~/PycharmProjects/python3projects/gtex_finemapping/ems_geuvadis_validation/train_a_predictor_for_a_chr_downstream.py gs://qingbowang/
#vmで:
gsutil cp gs://qingbowang/train_a_predictor_for_a_chr_downstream.py ./


#with multiple vms:

python3 train_a_predictor_for_a_chr_downstream.py chr1 ;
python3 train_a_predictor_for_a_chr_downstream.py chr2 ;
python3 train_a_predictor_for_a_chr_downstream.py chr3 ;
python3 train_a_predictor_for_a_chr_downstream.py chr4 #vm
python3 train_a_predictor_for_a_chr_downstream.py chr5 ;
python3 train_a_predictor_for_a_chr_downstream.py chr6 ;
python3 train_a_predictor_for_a_chr_downstream.py chr7 ;
python3 train_a_predictor_for_a_chr_downstream.py chr8 #vm2
python3 train_a_predictor_for_a_chr_downstream.py chr9 ;
python3 train_a_predictor_for_a_chr_downstream.py chr10 ;
python3 train_a_predictor_for_a_chr_downstream.py chr11 ;
python3 train_a_predictor_for_a_chr_downstream.py chr12 #vm3
python3 train_a_predictor_for_a_chr_downstream.py chr13 ;
python3 train_a_predictor_for_a_chr_downstream.py chr14 ;
python3 train_a_predictor_for_a_chr_downstream.py chr15 ;
python3 train_a_predictor_for_a_chr_downstream.py chr16 #vm4
python3 train_a_predictor_for_a_chr_downstream.py chr17 ;
python3 train_a_predictor_for_a_chr_downstream.py chr18 ;
python3 train_a_predictor_for_a_chr_downstream.py chr19 ;
python3 train_a_predictor_for_a_chr_downstream.py chr20 #vm5
python3 train_a_predictor_for_a_chr_downstream.py chr21 ;
python3 train_a_predictor_for_a_chr_downstream.py chr22 ;
python3 train_a_predictor_for_a_chr_downstream.py chrX #vm6
