#!/usr/bin/env bash

#the lines with "##" are done in local, whereas the others are done in a VM environment in google cloud

##cp the scripts
##gsutil cp ~/PycharmProjects/python3projects/gtex_finemapping/ems_pipe_test_20200125/train_a_predictor.py gs://qingbowang/
##gsutil cp ~/PycharmProjects/python3projects/gtex_finemapping/ems_pipe_test_20200125/train_a_predictor_for_a_tissue.sh gs://qingbowang/

##set up the VM
##gcloud beta compute instances create "rf-vm" --machine-type=n1-highcpu-96 --scopes=storage-rw
##gcloud beta compute --project "encode-uk-biobank" ssh --zone "us-central1-b" "rf-vm"


#load the packages to use
sudo apt-get install bzip2 libxml2-dev
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
pip install --user sklearn
pip install --user pandas
pip install --user matplotlib

#copy the script from cloud bucket
gsutil cp gs://qingbowang/train_a_predictor.py ./
gsutil cp gs://qingbowang/train_a_predictor_for_a_tissue.sh ./

#copy the files from cloud bucket
gsutil cp gs://qingbowang/ems_v1_test/bas_feature_top100_touse.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/roadmapannot_feat_to_use_all.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/n_pos_and_neg_per_tissue.tsv ./

#train a predictor per tissue
bash ./train_a_predictor_for_a_tissue.sh Whole_Blood 
bash ./train_a_predictor_for_a_tissue.sh Muscle_Skeletal 
bash ./train_a_predictor_for_a_tissue.sh Prostate 
bash ./train_a_predictor_for_a_tissue.sh Spleen
bash ./train_a_predictor_for_a_tissue.sh Skin_Sun_Exposed_Lower_leg 
bash ./train_a_predictor_for_a_tissue.sh Artery_Coronary 
bash ./train_a_predictor_for_a_tissue.sh Esophagus_Muscularis 
bash ./train_a_predictor_for_a_tissue.sh Esophagus_Gastroesophageal_Junction 
bash ./train_a_predictor_for_a_tissue.sh Artery_Tibial 
bash ./train_a_predictor_for_a_tissue.sh Heart_Atrial_Appendage 
bash ./train_a_predictor_for_a_tissue.sh Nerve_Tibial 
bash ./train_a_predictor_for_a_tissue.sh Heart_Left_Ventricle
bash ./train_a_predictor_for_a_tissue.sh Adrenal_Gland
bash ./train_a_predictor_for_a_tissue.sh Liver 
bash ./train_a_predictor_for_a_tissue.sh Adipose_Visceral_Omentum
bash ./train_a_predictor_for_a_tissue.sh Pancreas 
bash ./train_a_predictor_for_a_tissue.sh Lung
bash ./train_a_predictor_for_a_tissue.sh Pituitary
bash ./train_a_predictor_for_a_tissue.sh Brain_Nucleus_accumbens_basal_ganglia
bash ./train_a_predictor_for_a_tissue.sh Colon_Transverse
bash ./train_a_predictor_for_a_tissue.sh Adipose_Subcutaneous
bash ./train_a_predictor_for_a_tissue.sh Esophagus_Mucosa
bash ./train_a_predictor_for_a_tissue.sh Brain_Cerebellum
bash ./train_a_predictor_for_a_tissue.sh Brain_Cortex
bash ./train_a_predictor_for_a_tissue.sh Thyroid
bash ./train_a_predictor_for_a_tissue.sh Stomach
bash ./train_a_predictor_for_a_tissue.sh Breast_Mammary_Tissue
bash ./train_a_predictor_for_a_tissue.sh Colon_Sigmoid
bash ./train_a_predictor_for_a_tissue.sh Skin_Not_Sun_Exposed_Suprapubic
bash ./train_a_predictor_for_a_tissue.sh Testis
bash ./train_a_predictor_for_a_tissue.sh Artery_Aorta
bash ./train_a_predictor_for_a_tissue.sh Brain_Amygdala
bash ./train_a_predictor_for_a_tissue.sh Brain_Anterior_cingulate_cortex_BA24
bash ./train_a_predictor_for_a_tissue.sh Brain_Caudate_basal_ganglia
bash ./train_a_predictor_for_a_tissue.sh Brain_Hippocampus
bash ./train_a_predictor_for_a_tissue.sh Brain_Hypothalamus
bash ./train_a_predictor_for_a_tissue.sh Brain_Putamen_basal_ganglia
bash ./train_a_predictor_for_a_tissue.sh Cells_EBV-transformed_lymphocytes
bash ./train_a_predictor_for_a_tissue.sh Kidney_Cortex
bash ./train_a_predictor_for_a_tissue.sh Minor_Salivary_Gland
bash ./train_a_predictor_for_a_tissue.sh Brain_Cerebellar_Hemisphere
bash ./train_a_predictor_for_a_tissue.sh Brain_Frontal_Cortex_BA9
bash ./train_a_predictor_for_a_tissue.sh Ovary
bash ./train_a_predictor_for_a_tissue.sh Brain_Spinal_cord_cervical_c-1
bash ./train_a_predictor_for_a_tissue.sh Brain_Substantia_nigra
bash ./train_a_predictor_for_a_tissue.sh Cells_Cultured_fibroblasts
bash ./train_a_predictor_for_a_tissue.sh Small_Intestine_Terminal_Ileum
bash ./train_a_predictor_for_a_tissue.sh Uterus
bash ./train_a_predictor_for_a_tissue.sh Vagina


##and delete the vm
##gcloud beta compute instances delete "rf-vm"