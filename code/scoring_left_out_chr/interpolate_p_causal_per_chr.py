# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

#do this in random vm
import subprocess
import pandas as pd
import numpy as np

tissue_name = "Whole_Blood"
chr_list = []
for i in range(1,22+1):
    chr_list.append("chr" + str(i))
chr_list.append("chrX")


r = pd.DataFrame(pd.Series(np.arange(0,1.000001, 0.001)))
r.columns = ["rf_score_bin"]
r.index = r.rf_score_bin
for chr in chr_list:
    fn = "gs://qingbowang/ems_v1_test/ems_p_causal_{0}_leftout_{1}.tsv".format(tissue_name, chr)
    subprocess.call(["gsutil", "cp", fn, "./"])
    pcausal = pd.read_csv("ems_p_causal_{0}_leftout_{1}.tsv".format(tissue_name, chr), sep="\t", index_col=0)
    j = r.join(pcausal).interpolate(limit_direction="both") #joined, and interpolated for nans
    j["confidence_gain_log10"] = np.log10(j.confidence_gain)#need to re-calculate the log10 for interpolated ones
    #check there's no NA
    print("num. NAs: {0}".format(j.isna().sum().sum()))
    j.to_csv("ems_p_causal_interpolated_{0}_leftout_{1}.tsv".format(tissue_name, chr), sep="\t")
    subprocess.call(["gsutil", "cp", "ems_p_causal_interpolated_{0}_leftout_{1}.tsv".format(tissue_name, chr), "gs://qingbowang/ems_v1_test/"])
