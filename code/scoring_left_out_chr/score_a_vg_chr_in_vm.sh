#!/usr/bin/env bash


gcloud beta compute instances create "vg-vm2" --machine-type=n1-highmem-96 --scopes=storage-rw --boot-disk-size "500GB"
gcloud beta compute --project "encode-uk-biobank" ssh --zone "us-central1-b" "vg-vm2"


#load python things
sudo apt-get install bzip2 libxml2-dev
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
pip install --user sklearn
pip install --user pandas
#and the non tissue specific but needed file
gsutil cp gs://qingbowang/ems_v1_test/rf_feat_to_use_ordered_per_leftout_chr.tsv ./

#and test these in python script

#this in local
#gsutil cp ~/PycharmProjects/python3projects/gtex_finemapping/ems_geuvadis_validation/score_a_vg_chr_in_vm.py gs://qingbowang/
#this in vm
gsutil cp gs://qingbowang/score_a_vg_chr_in_vm.py ./



#以下はfor each chrに書き換える. cpは一回でOK


#this in vm6
for chr in chr1 chr2 chr3 chr4
do
    gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/vgpair_and_funcannot_"$chr".tsv.bgz ./
    
    for h in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9ac H3K9me3
    do
        gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/roadmap_vartouse_"$h"_"$chr".tsv.bgz ./
    done
    
    for chk in 0to499 500to999 1000to1499 1500to1999 2000to2499 2500to2999 3000to3499 3500to3999 4000to4499 4500to4999 5000to5311
    do
        gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/sad_vartouse_col"$chk"_"$chr".tsv.bgz ./
    done
    #also the predictor:
    gsutil cp gs://qingbowang/ems_v1_test/ems_rfmodel_Whole_Blood_leftout_"$chr".sav ./
    #and score
    python3 score_a_vg_chr_in_vm.py "$chr"
done

#this in vm5
for chr in chr5 chr6 chr7 chr8
do
    gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/vgpair_and_funcannot_"$chr".tsv.bgz ./

    for h in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9ac H3K9me3
    do
        gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/roadmap_vartouse_"$h"_"$chr".tsv.bgz ./
    done

    for chk in 0to499 500to999 1000to1499 1500to1999 2000to2499 2500to2999 3000to3499 3500to3999 4000to4499 4500to4999 5000to5311
    do
        gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/sad_vartouse_col"$chk"_"$chr".tsv.bgz ./
    done
    #also the predictor:
    gsutil cp gs://qingbowang/ems_v1_test/ems_rfmodel_Whole_Blood_leftout_"$chr".sav ./
    #and score
    python3 score_a_vg_chr_in_vm.py "$chr"
done

#vm4
for chr in chr9 chr10 chr11 chr12 chr13
do
    gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/vgpair_and_funcannot_"$chr".tsv.bgz ./

    for h in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9ac H3K9me3
    do
        gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/roadmap_vartouse_"$h"_"$chr".tsv.bgz ./
    done

    for chk in 0to499 500to999 1000to1499 1500to1999 2000to2499 2500to2999 3000to3499 3500to3999 4000to4499 4500to4999 5000to5311
    do
        gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/sad_vartouse_col"$chk"_"$chr".tsv.bgz ./
    done
    #also the predictor:
    gsutil cp gs://qingbowang/ems_v1_test/ems_rfmodel_Whole_Blood_leftout_"$chr".sav ./
    #and score
    python3 score_a_vg_chr_in_vm.py "$chr"
done

#vm3
for chr in chr14 chr15 chr16 chr17 chr18
do
    gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/vgpair_and_funcannot_"$chr".tsv.bgz ./

    for h in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9ac H3K9me3
    do
        gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/roadmap_vartouse_"$h"_"$chr".tsv.bgz ./
    done

    for chk in 0to499 500to999 1000to1499 1500to1999 2000to2499 2500to2999 3000to3499 3500to3999 4000to4499 4500to4999 5000to5311
    do
        gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/sad_vartouse_col"$chk"_"$chr".tsv.bgz ./
    done
    #also the predictor:
    gsutil cp gs://qingbowang/ems_v1_test/ems_rfmodel_Whole_Blood_leftout_"$chr".sav ./
    #and score
    python3 score_a_vg_chr_in_vm.py "$chr"
done

#vm2
for chr in chr19 chr20 chr21 chr22 chrX
do
    gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/vgpair_and_funcannot_"$chr".tsv.bgz ./

    for h in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9ac H3K9me3
    do
        gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/roadmap_vartouse_"$h"_"$chr".tsv.bgz ./
    done

    for chk in 0to499 500to999 1000to1499 1500to1999 2000to2499 2500to2999 3000to3499 3500to3999 4000to4499 4500to4999 5000to5311
    do
        gsutil cp gs://qingbowang/ems_v1_test/bgz_chunks/sad_vartouse_col"$chk"_"$chr".tsv.bgz ./
    done
    #also the predictor:
    gsutil cp gs://qingbowang/ems_v1_test/ems_rfmodel_Whole_Blood_leftout_"$chr".sav ./
    #and score
    python3 score_a_vg_chr_in_vm.py "$chr"
done


#delete the vm
gcloud beta compute instances delete "vg-vm6" #and so on