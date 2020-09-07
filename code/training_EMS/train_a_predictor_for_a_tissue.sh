#!/usr/bin/env bash


tissue_name=$1

#copy the files from the bucket
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_annotated.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_annotated.tsv ./

gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_roadmapannot_H3K27ac.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_roadmapannot_H3K27me3.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_roadmapannot_H3K36me3.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_roadmapannot_H3K4me1.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_roadmapannot_H3K4me3.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_roadmapannot_H3K9ac.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_positive_training_vg_roadmapannot_H3K9me3.tsv ./

gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_roadmapannot_H3K27ac.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_roadmapannot_H3K27me3.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_roadmapannot_H3K36me3.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_roadmapannot_H3K4me1.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_roadmapannot_H3K4me3.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_roadmapannot_H3K9ac.tsv ./
gsutil cp gs://qingbowang/ems_v1_test/"$tissue_name"_negative_training_vg_roadmapannot_H3K9me3.tsv ./

gsutil cp gs://qingbowang/ems_v1_test/n_pos_and_neg_per_tissue.tsv ./ #when we also do the save the predictor part

#parameter search:
python3 train_a_predictor.py $tissue_name

#cp the product file
gsutil cp ./grid_search_result_"$tissue_name".tsv gs://qingbowang/ems_v1_test/grid_search_result_"$tissue_name".tsv
gsutil cp ./random_search_result_"$tissue_name".tsv gs://qingbowang/ems_v1_test/random_search_result_"$tissue_name".tsv

#rm the file
rm ./"$tissue_name"_positive_training_vg_annotated.tsv
rm ./"$tissue_name"_negative_training_vg_annotated.tsv

rm ./"$tissue_name"_positive_training_vg_roadmapannot_H3K27ac.tsv
rm ./"$tissue_name"_positive_training_vg_roadmapannot_H3K27me3.tsv
rm ./"$tissue_name"_positive_training_vg_roadmapannot_H3K36me3.tsv
rm ./"$tissue_name"_positive_training_vg_roadmapannot_H3K4me1.tsv
rm ./"$tissue_name"_positive_training_vg_roadmapannot_H3K4me3.tsv
rm ./"$tissue_name"_positive_training_vg_roadmapannot_H3K9ac.tsv
rm ./"$tissue_name"_positive_training_vg_roadmapannot_H3K9me3.tsv
rm ./"$tissue_name"_negative_training_vg_roadmapannot_H3K27ac.tsv
rm ./"$tissue_name"_negative_training_vg_roadmapannot_H3K27me3.tsv
rm ./"$tissue_name"_negative_training_vg_roadmapannot_H3K36me3.tsv
rm ./"$tissue_name"_negative_training_vg_roadmapannot_H3K4me1.tsv
rm ./"$tissue_name"_negative_training_vg_roadmapannot_H3K4me3.tsv
rm ./"$tissue_name"_negative_training_vg_roadmapannot_H3K9ac.tsv
rm ./"$tissue_name"_negative_training_vg_roadmapannot_H3K9me3.tsv
