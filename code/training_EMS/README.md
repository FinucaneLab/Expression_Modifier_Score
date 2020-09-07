This directory contains codes used to train a predictor to build EMS. 

`prep_training_data.sh` was used to annotate the GTEx variant-genes (vg) with posterior probability and functional features.
The code calls
- `bin_annotate_gtex_var.py` which annotates GTEx vg with baseline annotations, 
- `write_gtex_assoc_into_ht.py` which writes out the GTEx vg in hail table, 
- `annotate_bin_and_pp.py` which annotates the vg with posterior probability,
- `export_pos_and_neg_for_training.py` which further annotates the positive/negative training data with ROADMAP and Basenji features, 
- `select_tissue_spec_bin_feat.py` which selects the ROADMAP functional features that will be used for the predictor, based on F1 measure etc, and 
- `select_bas_feat.py` which selects Basenji features that will be used for the predictor, based on marginal ROC.

`select_tissue_spec_bin_feat_downstream.py` and `select_bas_feat_downstream.py` were used to summarized the results across tissue types.

`train_a_predictor_for_a_tissue.sh` (which calls `train_a_predictor.py`) was used to optimize the parameters in a random forest (RF) predictor that is used for building EMS by performing random and grid parameter search.

`train_predictors_alltissues.sh` was used to run the code above for all 49 GTEx tissues.

`train_predictors_alltissues_downstream.sh`, (which calls `train_a_predictor_downstream.py` for all tissues) was used to train the predictor with optimized parameter.

Finally, `output_features_to_use_ordered.py` and `interpolate_p_causal.py` were used to format the files so that they can be streamlined for the downstream processes (e.g. scoring all the GTEx variants using the RF)
