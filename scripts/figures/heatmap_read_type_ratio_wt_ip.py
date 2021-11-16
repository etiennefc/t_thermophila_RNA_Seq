#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

wt_df = pd.read_csv(snakemake.input.WT_df, sep='\t')
ip_df = pd.read_csv(snakemake.input.MLP1_IP_df, sep='\t')

wt_df = wt_df.drop(columns=['other_reads_norm']).set_index('isotype_anticodon')
ip_df = ip_df.drop(columns=['other_reads_norm']).set_index('isotype_anticodon')

# Multiply all values by 10^6 to get a better representation on the heatmap
wt_df, ip_df = wt_df * 1000000, ip_df * 1000000

# Change 0 to 0.001 so that log10 can be applied on df
wt_df = wt_df.replace(0, 0.001)
ip_df = ip_df.replace(0, 0.001)

# Divide value of ip_df by wt_df
div_df = ip_df.div(wt_df)

# Create simpler column names and get log10 of all values
cols = []
for col in list(div_df.columns):
    simple_col = col.split('_')[1]
    cols.append(simple_col)
    div_df[col+'_log10'] = np.log10(div_df[col])

div_df = div_df.filter(like='_log10', axis=1)

# Create heatmap (not clustered for isotypes nor read_type)
ft.heatmap_fixed(div_df, 'plasma', '\nlog10(average IP/WT RPM ratio)', cols,
                snakemake.output.heatmap, col_cluster=False, row_cluster=False)
