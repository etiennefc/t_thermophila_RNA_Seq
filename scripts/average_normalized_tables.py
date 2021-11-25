#!/usr/bin/python3
import pandas as pd

# Load all dfs per condition
wt, ko, ip = [], [], []
for df_path in snakemake.input.normalized_dfs:
    if 'R1' in df_path:  # select only reads R1 (not R2)
        if 'WT' in df_path:
            df = pd.read_csv(df_path, sep='\t')
            df = df.fillna(0)
            wt.append(df)
        elif 'KO' in df_path:
            df = pd.read_csv(df_path, sep='\t')
            df = df.fillna(0)
            ko.append(df)
        else:
            df = pd.read_csv(df_path, sep='\t')
            df = df.fillna(0)
            ip.append(df)

# Concat dfs per condition
wt_concat, ko_concat, ip_concat = pd.concat(wt), pd.concat(ko), pd.concat(ip)

# Get an average value per condition
wt_avg = wt_concat.groupby('isotype_anticodon').mean()
ko_avg = ko_concat.groupby('isotype_anticodon').mean()
ip_avg = ip_concat.groupby('isotype_anticodon').mean()

# Select normalized columns only and isotype_anticodon
wt_avg = wt_avg.filter(like='_norm', axis=1).reset_index()
ko_avg = ko_avg.filter(like='_norm', axis=1).reset_index()
ip_avg = ip_avg.filter(like='_norm', axis=1).reset_index()

wt_avg.to_csv(snakemake.output.WT_df, index=False, sep='\t')
ko_avg.to_csv(snakemake.output.KO_df, index=False, sep='\t')
ip_avg.to_csv(snakemake.output.MLP1_IP_df, index=False, sep='\t')


