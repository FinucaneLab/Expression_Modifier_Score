# Expression Modifier Score (EMS)
Welcome to the github page for the ["Predicting the probability of cis-regulatory effects in the human genome with supervised learning from 14,807 fine-mapped eQTLs and 6,121 features"](link) manuscript, which introduces the Expression Modifier Score (EMS). 

The `tutorials` directory contains all the information necessary to annotate EMS on your own dataset.
It also allows user to reproduce the main figure and most of the supplementary figure in the manuscript. 
Specifically, the tutorials directory consists of xx Jupyter notebooks:
1. `annotate_EMS.ipynb` explains how to annotate EMS on your dataset.
2. xx
3. xx


Which figure and table in the preprint/paper is generated in which notebook is listed below:
(those with () does are related but not exactly used)

| notebook  | figure  |   supplementary file  |
|---|---|---|
|annotate_EMS.ipynb   | xx  | xx  |
|xx.ipynb   | xx  |   |
|xx.ipynb  | xx  | xx  |
|xx.ipynb  | xx  | xx  |
|xx.ipynb  | xx  |  |
|xx.ipynb  |   |  |

Since most of the analysis was performed in Hail, we recommend users who are not familiat with Hail to visit the [Hail tutorial page](https://hail.is/docs/0.2/tutorials-landing.html).
(The analysis was performed using Hail version 0.2.11, and we recommend downloading this specific version of hail to perform MNV analysis using Hail, e.g. with command `pip install hail==0.2.11`.)


All the scripts used in the EMS paper are stored in the `code` directory. 
However, note that due to some of the genomics data being not publicly available, 
most of the scripts cannot be simply run in your local.

`util` contains the functions used in the analysis.
