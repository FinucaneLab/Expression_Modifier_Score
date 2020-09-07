#!/usr/bin/env bash


tissue_name=Whole_Blood
leftout_chr=$1

#cp the files
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


#run the python script for param optimization (faster if parralel, but fine for now...)
python3 train_a_predictor_for_a_chr.py $tissue_name $leftout_chr

#cp the product file
gsutil cp ./grid_search_result_"$tissue_name"_"$leftout_chr"_left_out.tsv gs://qingbowang/ems_v1_test/grid_search_result_"$leftout_chr".tsv
gsutil cp ./random_search_result_"$tissue_name"_"$leftout_chr"_left_out.tsv gs://qingbowang/ems_v1_test/random_search_result_"$leftout_chr".tsv



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
