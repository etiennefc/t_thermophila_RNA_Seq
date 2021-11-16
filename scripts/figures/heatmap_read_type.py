#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

wt_df = pd.read_csv(snakemake.input.WT_df, sep='\t')
ko_df = pd.read_csv(snakemake.input.KO_df, sep='\t')
ip_df = pd.read_csv(snakemake.input.MLP1_IP_df, sep='\t')

wt_df = wt_df.drop(columns=['other_reads_norm']).set_index('isotype_anticodon')
ko_df = ko_df.drop(columns=['other_reads_norm']).set_index('isotype_anticodon')
ip_df = ip_df.drop(columns=['other_reads_norm']).set_index('isotype_anticodon')

# Multiply all values by 10^6 to get a better representation on the heatmap
wt_df, ko_df, ip_df = wt_df * 1000000, ko_df * 1000000, ip_df * 1000000

# Change 0 to 0.001 so that log10 can be applied on df
wt_df = wt_df.replace(0, 0.001)
ko_df = ko_df.replace(0, 0.001)
ip_df = ip_df.replace(0, 0.001)

# Create simpler column names and get log10 of all values
cols = []
for col in list(wt_df.columns):
    simple_col = col.split('_')[1]
    cols.append(simple_col)
    wt_df[col+'_log10'] = np.log10(wt_df[col])
    ko_df[col+'_log10'] = np.log10(ko_df[col])
    ip_df[col+'_log10'] = np.log10(ip_df[col])

wt_df = wt_df.filter(like='_log10', axis=1)
ko_df = ko_df.filter(like='_log10', axis=1)
ip_df = ip_df.filter(like='_log10', axis=1)

# Create heatmap (clustered for isotypes, but not read  type)
ft.heatmap_fixed(wt_df, 'plasma', '\nlog10(average RPM)', cols,
                snakemake.output.heatmap_wt, col_cluster=False, row_cluster=False, method='weighted')
ft.heatmap_fixed(ko_df, 'plasma', '\nlog10(average RPM)', cols,
                snakemake.output.heatmap_ko, col_cluster=False, row_cluster=False, method='weighted')
ft.heatmap_fixed(ip_df, 'plasma', '\nlog10(average RPM)', cols,
                snakemake.output.heatmap_ip, col_cluster=False, row_cluster=False, method='weighted')
