#!/usr/bin/env bash


#the lines with "##" are done in local, whereas the others are done in a VM environment in google cloud

##cp the scripts
##gsutil cp ~/PycharmProjects/python3projects/gtex_finemapping/ems_pipe_test_20200125/train_a_predictor_downstream.py gs://qingbowang/

##set up the VM
##gcloud beta compute instances create "rf-vm" --machine-type=n1-highcpu-96 --scopes=storage-rw
##gcloud beta compute --project "encode-uk-biobank" ssh --zone "us-central1-b" "rf-vm"




gsutil cp gs://qingbowang/train_a_predictor_downstream.py ./

python3 ./train_a_predictor_downstream.py Whole_Blood 
python3 ./train_a_predictor_downstream.py Muscle_Skeletal 
python3 ./train_a_predictor_downstream.py Prostate 
python3 ./train_a_predictor_downstream.py Spleen
python3 ./train_a_predictor_downstream.py Skin_Sun_Exposed_Lower_leg 
python3 ./train_a_predictor_downstream.py Artery_Coronary 
python3 ./train_a_predictor_downstream.py Esophagus_Muscularis 
python3 ./train_a_predictor_downstream.py Esophagus_Gastroesophageal_Junction 
python3 ./train_a_predictor_downstream.py Artery_Tibial 
python3 ./train_a_predictor_downstream.py Heart_Atrial_Appendage 
python3 ./train_a_predictor_downstream.py Nerve_Tibial 
python3 ./train_a_predictor_downstream.py Heart_Left_Ventricle
python3 ./train_a_predictor_downstream.py Adrenal_Gland
python3 ./train_a_predictor_downstream.py Liver 
python3 ./train_a_predictor_downstream.py Adipose_Visceral_Omentum
python3 ./train_a_predictor_downstream.py Pancreas 
python3 ./train_a_predictor_downstream.py Lung
python3 ./train_a_predictor_downstream.py Pituitary
python3 ./train_a_predictor_downstream.py Brain_Nucleus_accumbens_basal_ganglia
python3 ./train_a_predictor_downstream.py Colon_Transverse
python3 ./train_a_predictor_downstream.py Adipose_Subcutaneous
python3 ./train_a_predictor_downstream.py Esophagus_Mucosa
python3 ./train_a_predictor_downstream.py Brain_Cerebellum
python3 ./train_a_predictor_downstream.py Brain_Cortex
python3 ./train_a_predictor_downstream.py Thyroid
python3 ./train_a_predictor_downstream.py Stomach
python3 ./train_a_predictor_downstream.py Breast_Mammary_Tissue
python3 ./train_a_predictor_downstream.py Colon_Sigmoid
python3 ./train_a_predictor_downstream.py Skin_Not_Sun_Exposed_Suprapubic
python3 ./train_a_predictor_downstream.py Testis
python3 ./train_a_predictor_downstream.py Artery_Aorta
python3 ./train_a_predictor_downstream.py Brain_Amygdala
python3 ./train_a_predictor_downstream.py Brain_Anterior_cingulate_cortex_BA24
python3 ./train_a_predictor_downstream.py Brain_Caudate_basal_ganglia
python3 ./train_a_predictor_downstream.py Brain_Hippocampus
python3 ./train_a_predictor_downstream.py Brain_Hypothalamus
python3 ./train_a_predictor_downstream.py Brain_Putamen_basal_ganglia
python3 ./train_a_predictor_downstream.py Cells_EBV-transformed_lymphocytes
python3 ./train_a_predictor_downstream.py Kidney_Cortex
python3 ./train_a_predictor_downstream.py Minor_Salivary_Gland
python3 ./train_a_predictor_downstream.py Brain_Cerebellar_Hemisphere
python3 ./train_a_predictor_downstream.py Brain_Frontal_Cortex_BA9
python3 ./train_a_predictor_downstream.py Ovary
python3 ./train_a_predictor_downstream.py Brain_Spinal_cord_cervical_c-1
python3 ./train_a_predictor_downstream.py Brain_Substantia_nigra
python3 ./train_a_predictor_downstream.py Cells_Cultured_fibroblasts
python3 ./train_a_predictor_downstream.py Small_Intestine_Terminal_Ileum
python3 ./train_a_predictor_downstream.py Uterus
python3 ./train_a_predictor_downstream.py Vagina


##and delete the vm
##gcloud beta compute instances delete "rf-vm"