This directory contains codes used to score GTEx variant-gene pairs (vg) with EMS.

- `score_gtex_vg_prep.sh` (which calls `score_getx_vg_prep.py`) was used to format the data.

- `score_gtex_vg_alltissues.sh` (which calls `score_a_vg_tissue_in_vm.py` for all the tissues) was used to score the GTEx vg, and 

- `score_gtex_vg_downstream.py` to format the output and save as hail table.

When evaluating the performance of EMS (as in the fig. 1 of the manuscript), each chromosome was scored using the predictor that was trained excluding the chromosome, to avoid overfitting. 
Those codes can be found in the `scoring_left_out_chr` directory.
