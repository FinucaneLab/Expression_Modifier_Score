#!/usr/bin/env bash

#the lines with "##" are done in local, whereas the others are done in a VM environment in google cloud
##gcloud beta compute instances create "vg-vm" --machine-type=n1-highmem-96 --scopes=storage-rw --boot-disk-size "500GB"
##gcloud beta compute --project "encode-uk-biobank" ssh --zone "us-central1-b" "vg-vm"
##gsutil cp ~/PycharmProjects/python3projects/gtex_finemapping/ems_pipe_test_20200125/score_a_vg_tissue_in_vm.py gs://qingbowang/


#load other necessary files/packages
sudo apt-get install bzip2 libxml2-dev
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
pip install --user sklearn
pip install --user pandas
gsutil cp gs://qingbowang/ems_v1_test/rf_feat_to_use_ordered.tsv ./
gsutil cp gs://qingbowang/score_a_vg_tissue_in_vm.py ./

#do the scoring for all tissues




for tissue_name in Whole_Blood Muscle_Skeletal Liver Brain_Cerebellum Prostate Spleen Skin_Sun_Exposed_Lower_leg Artery_Coronary Esophagus_Muscularis Esophagus_Gastroesophageal_Junction Artery_Tibial Heart_Atrial_Appendage Nerve_Tibial Heart_Left_Ventricle Adrenal_Gland Adipose_Visceral_Omentum Pancreas Lung Pituitary Brain_Nucleus_accumbens_basal_ganglia Colon_Transverse Adipose_Subcutaneous Esophagus_Mucosa Brain_Cortex Thyroid Stomach Breast_Mammary_Tissue Colon_Sigmoid Skin_Not_Sun_Exposed_Suprapubic Testis Artery_Aorta Brain_Amygdala Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus Brain_Putamen_basal_ganglia Brain_Spinal_cord_cervical_c-1 Brain_Substantia_nigra Cells_Cultured_fibroblasts Cells_EBV-transformed_lymphocytes Kidney_Cortex Minor_Salivary_Gland Ovary Small_Intestine_Terminal_Ileum Uterus Vagina
do
    #load necessary files
    gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/vgpair_and_funcannot_"$tissue_name".tsv.bgz ./
    
    for h in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9ac H3K9me3
    do
        gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/roadmap_vartouse_"$h"_"$tissue_name".tsv.bgz ./
    done
    
    for chk in 0to499 500to999 1000to1499 1500to1999 2000to2499 2500to2999 3000to3499 3500to3999 4000to4499 4500to4999 5000to5311
    do
        gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/sad_vartouse_col"$chk"_"$tissue_name".tsv.bgz ./
    done
    #also load the predictor:
    gsutil cp gs://qingbowang/ems_v1_test/ems_rfmodel_"$tissue_name".sav ./
    #and score
    python3 score_a_vg_tissue_in_vm.py "$tissue_name"
done

#delete the vm
##gcloud beta compute instances delete "vg-vm"