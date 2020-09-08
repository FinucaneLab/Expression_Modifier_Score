#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 12:14:34 2020

@author: shuvomsadhuka
"""
from sklearn import metrics
import matplotlib.pyplot as plt
import math
import numpy as np
import os
import pandas as pd

plt.rcParams.update({'font.size': 16})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
pd.set_option('display.max_columns', None)


def DeepSea(df):
    e_vals = df.iloc[:, 8:]
    #e_vals['min'] = e_vals.min(axis=1, skipna=False, numeric_only=True)
    e_vals = e_vals.fillna(0)
    e_vals = e_vals.mask(e_vals==0).fillna(e_vals.min(axis=1))
    e_vals = -np.log10(e_vals)
    e_vals['mean e_val'] = e_vals.mean(axis=1)
    e_vals = e_vals.dropna()
    #e_vals['POS'], e_vals['REF'], e_vals['ALT'] = df['pos'], df['ref'], df['alt']
    return(e_vals)


def ems(df, plt_type):
    unravel_snp = df['vg'].str.split('_', expand = True)
    df['CHROM'], df['POS'], df['REF'], df['ALT'] = unravel_snp[0], unravel_snp[1], unravel_snp[2], unravel_snp[3]
    df_max = df.sort_values('ems', 
                            ascending=False).drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT']).sort_index()
    num_pos = len(df[df['label'] == 1])
    num_neg = len(df[df['label'] == 0])
    df['sample_weight'] = np.where(df['label']==1, 1, num_pos/num_neg)
    
    num_pos = len(df_max[df_max['label'] == 1])
    num_neg = len(df_max[df_max['label'] == 0])
    df_max['sample_weight'] = np.where(df_max['label']==1, 1, num_pos/num_neg)
    
    if plt_type == 'prc':
        precision, recall, thresholds = metrics.precision_recall_curve(df_max['label'], df_max['ems'])
        auc = metrics.average_precision_score(df_max['label'], df_max['ems'])
        plt.plot(recall, precision, c='tab:green', label = 'EMS, AUC = %0.3f' % auc)
        
        precision, recall, thresholds = metrics.precision_recall_curve(df['label'], df['ems'])
        auc = metrics.average_precision_score(df['label'], df['ems'])
        #plt.plot(recall, precision, c='tab:green', label = 'EMS, AUC = %0.3f' % auc)
        
    if plt_type == 'roc':
        fpr, tpr, threshold = metrics.roc_curve(df_max['label'], df_max['ems'])
        roc_auc = metrics.auc(fpr, tpr)
        #plt.plot(fpr, tpr, c='tab:green', label = 'EMS, AUC = %0.3f' % roc_auc)
        
        fpr, tpr, threshold = metrics.roc_curve(df['label'], df['ems'], sample_weight=df['sample_weight'])
        roc_auc = metrics.auc(fpr, tpr)
        plt.plot(fpr, tpr, c='tab:green', label = 'EMS, AUC = %0.3f' % roc_auc)
        

def make_prc(df_pos, df_neg, version):
    df_pos = pd.read_csv(df_pos, index_col=[0])
    df_pos['label'] = [1 for i in range(len(df_pos))]
    #print(df_pos[df_pos.duplicated(subset=['v'])])
    
    df_neg = pd.read_csv(df_neg, index_col=[0])
    df_neg['label'] = [0 for i in range(len(df_neg))]
    df = pd.concat([df_pos, df_neg])
    df.rename(columns={'vg_x': 'vg'}, inplace=True)
    
    plt.title('Per Variant Precision Recall Curve')
    if version == 'cadd':
        df = df.dropna(subset=[version])
        unravel_snp = df['vg'].str.split('_', expand = True)
        df['CHROM'], df['POS'], df['REF'], df['ALT'] = unravel_snp[0], unravel_snp[1], unravel_snp[2], unravel_snp[3]
        df = df.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT']).sort_index()
        num_pos = len(df[df['label'] == 1])
        num_neg = len(df[df['label'] == 0])
        df['sample_weight'] = np.where(df['label']==1, 1, num_pos/num_neg)
        precision, recall, thresholds = metrics.precision_recall_curve(df['label'], df[version])
        auc = metrics.average_precision_score(df['label'], df[version])
        plt.plot(recall, precision, c='tab:pink', label = 'CADD, AUC = %0.3f' % auc)
        
    if version == 'gerp':
        df = df.dropna(subset=[version])
        unravel_snp = df['vg'].str.split('_', expand = True)
        df['CHROM'], df['POS'], df['REF'], df['ALT'] = unravel_snp[0], unravel_snp[1], unravel_snp[2], unravel_snp[3]
        df = df.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT']).sort_index()
        df[version] = abs(df[version])
        num_pos = len(df[df['label'] == 1])
        num_neg = len(df[df['label'] == 0])
        df['sample_weight'] = np.where(df['label']==1, 1, num_pos/num_neg)
        precision, recall, thresholds = metrics.precision_recall_curve(df['label'], abs(df[version]))
        auc = metrics.average_precision_score(df['label'], df[version])
        plt.plot(recall, precision, c='tab:cyan', label = 'GERP, AUC = %0.3f' % auc)
        
    if version == 'ems':
        df = df.dropna(subset=[version])
        ems(df, 'prc')
        
    if version == 'tss':
        df = df[df['tss_distance'] != 0]
        df['tss'] = 1/abs(df['tss_distance'])
        df = df.dropna(subset=[version])
        unravel_snp = df['vg'].str.split('_', expand = True)
        df['CHROM'], df['POS'], df['REF'], df['ALT'] = unravel_snp[0], unravel_snp[1], unravel_snp[2], unravel_snp[3]
        df_max = df.sort_values('tss', 
                            ascending=False).drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT']).sort_index()
        num_pos = len(df_max[df_max['label'] == 1])
        num_neg = len(df_max[df_max['label'] == 0])
        df['sample_weight'] = np.where(df['label']==1, 1, num_pos/num_neg)
        precision, recall, thresholds = metrics.precision_recall_curve(df_max['label'], abs(df_max[version]))
        auc = metrics.average_precision_score(df_max['label'], df_max[version])
        plt.plot(recall, precision, c='tab:orange', label = 'TSS, AUC = %0.3f' % auc)
        
    plt.legend(loc = 'lower right')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.xlabel('Recall', fontsize=20)
    plt.ylabel('Precision',  fontsize=20)
    

plt.figure(figsize=(7.5,5))

make_prc('Annotated/12k_pos.tsv', 'Annotated/12k_neg.tsv', 'ems')
make_prc('Annotated/12k_pos.tsv', 'Annotated/12k_neg.tsv', 'tss')

# deepsea plot
deepsea_pos = pd.read_csv('Results/deepsea_pos_12k/deepsea_pos_12k.tsv',
                      sep='\t')

deepsea_neg = pd.read_csv('Results/deepsea_neg_12k/deepsea_neg_12k.tsv', sep='\t')

deepsea_pos = DeepSea(deepsea_pos)
deepsea_pos['label'] = [1 for i in range(len(deepsea_pos))]
deepsea_neg = DeepSea(deepsea_neg)
deepsea_neg['label'] = [0 for i in range(len(deepsea_neg))]

df = pd.concat([deepsea_pos, deepsea_neg])
df = df.drop_duplicates(subset=['8988T|DNase|None', 'AoSMC|DNase|None'])
num_pos = len(df[df['label'] == 1])
num_neg = len(df[df['label'] == 0])
df['sample_weight'] = np.where(df['label']==1, 1, num_pos/num_neg)


fpr, tpr, threshold = metrics.roc_curve(df['label'], df['mean e_val'])
roc_auc = metrics.auc(fpr, tpr)
#plt.plot(fpr, tpr, c='tab:blue', label = 'DeepSEA, AUC = %0.3f' % roc_auc)

precision, recall, thresholds = metrics.precision_recall_curve(df['label'], df['mean e_val'])
auc = metrics.average_precision_score(df['label'], df['mean e_val'])
plt.plot(recall, precision, c='tab:blue', label = 'DeepSEA, AUC = %0.3f' % auc)


make_prc('Annotated/12k_pos.tsv', 'Annotated/12k_neg.tsv', 'cadd')

## NCER
ncer_pos = pd.read_csv('Annotated/12k_pos_ncer.tsv', sep='\t')
ncer_pos['label'] = [1 for i in range(len(ncer_pos))]
ncer_neg = pd.read_csv('Annotated/12k_neg_ncer.tsv', sep='\t')
ncer_neg['label'] = [0 for i in range(len(ncer_neg))]
ncer_df = pd.concat([ncer_pos, ncer_neg])
ncer_df = ncer_df.dropna()
ncer_df = ncer_df.drop_duplicates(subset=['variant_id'])
fpr, tpr, threshold = metrics.roc_curve(ncer_df['label'], pd.to_numeric(ncer_df['ncer']))
roc_auc = metrics.auc(fpr, tpr)
#plt.plot(fpr, tpr, c='tab:olive', label = 'ncER, AUC = %0.3f' % roc_auc)

precision, recall, thresholds = metrics.precision_recall_curve(ncer_df['label'], 
                                                              pd.to_numeric(ncer_df['ncer']))
auc = metrics.average_precision_score(ncer_df['label'], pd.to_numeric(ncer_df['ncer']))
plt.plot(recall, precision, c='tab:olive', label = 'ncER, AUC = %0.3f' % auc)

make_prc('Annotated/12k_pos.tsv', 'Annotated/12k_neg.tsv', 'gerp')


# fathmm plot
fathmm_neg = pd.read_csv('Results/final_results/fathmm_12k_neg.txt', sep='\t')
fathmm_neg['label'] = [0 for i in range(len(fathmm_neg))]

fathmm_pos = pd.read_csv('Results/final_results/fathmm_12k_pos.txt', sep='\t')
fathmm_pos['label'] = [1 for i in range(len(fathmm_pos))]

fathmm_df = pd.concat([fathmm_pos, fathmm_neg])
fathmm_df = fathmm_df[pd.to_numeric(fathmm_df['Non-Coding Score'], errors='coerce').notnull()]
fathmm_df = fathmm_df.drop_duplicates(subset=['Position', 'Ref. Base', 'Mutant Base'])
num_pos = len(fathmm_df[fathmm_df['label'] == 1])
num_neg = len(fathmm_df[fathmm_df['label'] == 0])
fathmm_df['sample_weight'] = np.where(fathmm_df['label']==1, 1, num_pos/num_neg)

fpr, tpr, threshold = metrics.roc_curve(fathmm_df['label'], pd.to_numeric(fathmm_df['Non-Coding Score']))
roc_auc = metrics.auc(fpr, tpr)
#plt.title('Receiver Operating Characteristic')
#plt.plot(fpr, tpr, c='tab:brown', label = 'Fathmm, AUC = %0.3f' % roc_auc)
precision, recall, thresholds = metrics.precision_recall_curve(fathmm_df['label'], 
                                                              pd.to_numeric(fathmm_df['Non-Coding Score']))
auc = metrics.average_precision_score(fathmm_df['label'], pd.to_numeric(fathmm_df['Non-Coding Score']))
plt.plot(recall, precision, c='tab:brown', label = 'Fathmm, AUC = %0.3f' % auc)






plt.legend(loc='upper right', prop={'size': 14})

file_loc = 'Plots/all_vg_prc_12k.png'
#plt.savefig(file_loc, dpi=400)