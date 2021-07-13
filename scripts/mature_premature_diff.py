#!/usr/bin/python3
import pandas as pd
import re


""" Get the difference between the number of reads associated to premature vs
    mature tRNAs per average sample."""

df_reads = pd.read_csv(snakemake.input.df, sep='\t')
sample_type = ['KO_', 'MLP1_IP_', 'WT_']

for i, sample in enumerate(sample_type):
    cols = list(df_reads.columns)

    #Select mature and preamture columns per sample and create average columns
    mature = [col for col in cols if "_mature" in col if sample in col]
    premature = [col for col in cols if "premature" in col if sample in col]

    df_reads[sample+'mature_avg'] = df_reads[mature].mean(axis=1)
    df_reads[sample+'premature_avg'] = df_reads[premature].mean(axis=1)

    # Create columns for the difference of number of reads (avg_premature - avg_mature per sample)
    df_reads[sample+'read_diff'] = df_reads[sample+'premature_avg'] - df_reads[sample+'mature_avg']


df_reads.to_csv(snakemake.output.avg_diff_df, index=False, sep='\t')
