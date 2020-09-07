#!/usr/bin/env bash

#cp the scripts
gsutil cp ~/PycharmProjects/python3projects/gtex_finemapping/ems_geuvadis_validation/train_a_predictor_for_a_chr.py gs://qingbowang/
gsutil cp ~/PycharmProjects/python3projects/gtex_finemapping/ems_geuvadis_validation/train_a_predictor_for_a_chr.sh gs://qingbowang/


gcloud beta compute instances create "rf-vm" --machine-type=n1-highcpu-96 --scopes=storage-rw
#gcloud beta compute instances create "rf-vm2" --machine-type=n2-highcpu-80 --scopes=storage-rw #not available again..
#gcloud beta compute instances create "rf-vm" --machine-type=c2-standard-60 --scopes=storage-rw #faster? no not available
gcloud beta compute --project "encode-uk-biobank" ssh --zone "us-central1-b" "rf-vm"
#and run this inside the ssh:

#load the package
sudo apt-get install bzip2 libxml2-dev
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
pip install --user sklearn
pip install --user pandas
#pip install --user pickle
pip install --user matplotlib

#copy the script from cloud
gsutil cp gs://qingbowang/train_a_predictor_for_a_chr.py ./
gsutil cp gs://qingbowang/train_a_predictor_for_a_chr.sh ./

#cp the files
gsutil cp gs://qingbowang/ems_v1_test/bas_feature_top100_touse_per_leftout_chr.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/roadmapannot_feat_to_use_per_leftout_chr.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/n_pos_and_neg_per_leftout_chr.tsv ./

#do them one by one (or maybe multiple clusters, in parallel)
bash ./train_a_predictor_for_a_chr.sh chr1
bash ./train_a_predictor_for_a_chr.sh chr2
bash ./train_a_predictor_for_a_chr.sh chr3
bash ./train_a_predictor_for_a_chr.sh chr4
bash ./train_a_predictor_for_a_chr.sh chr5 #1-5: in vm
bash ./train_a_predictor_for_a_chr.sh chr6 ;
bash ./train_a_predictor_for_a_chr.sh chr7 ;
bash ./train_a_predictor_for_a_chr.sh chr8 ;
bash ./train_a_predictor_for_a_chr.sh chr9 #6-9 in vm2
bash ./train_a_predictor_for_a_chr.sh chr10 ;
bash ./train_a_predictor_for_a_chr.sh chr11 ;
bash ./train_a_predictor_for_a_chr.sh chr12 ;
bash ./train_a_predictor_for_a_chr.sh chr13 #10-13 in vm3
bash ./train_a_predictor_for_a_chr.sh chr14 ;
bash ./train_a_predictor_for_a_chr.sh chr15 ;
bash ./train_a_predictor_for_a_chr.sh chr16 ;
bash ./train_a_predictor_for_a_chr.sh chr17 #14-17 in vm4
bash ./train_a_predictor_for_a_chr.sh chr18 ;
bash ./train_a_predictor_for_a_chr.sh chr19 ;
bash ./train_a_predictor_for_a_chr.sh chr20 #18-20 in vm5
bash ./train_a_predictor_for_a_chr.sh chr21 ;
bash ./train_a_predictor_for_a_chr.sh chr22 ;
bash ./train_a_predictor_for_a_chr.sh chrX #21,22,X in vm6


#and delete
gcloud beta compute instances delete "rf-vm"