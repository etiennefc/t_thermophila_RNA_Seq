#!/usr/bin/python3
import pandas as pd
import os

""" Generate a dataframe of gene biotype per gene_id from a gtf and merge that
    column to an abundance df. Add an average TPM column per sample type."""

cols = ['chr', 'source', 'feature', 'start', 'end', 'score1', 'strand', 'score2', 'features']
gtf = pd.read_csv(snakemake.input.gtf, sep='\t', names=cols)
gtf = gtf[gtf['feature'] == 'gene']
tpm_df = pd.read_csv(snakemake.input.tpm_df, sep='\t')

# Get gene biotype and note (gene function) for all genes in gtf from the features column
gtf_df = gtf[['features']]
gtf_df = gtf_df['features'].str.split(';', expand=True)
gtf_df.columns = ['gene_id', 'gene_name', 'gene_version', 'gene_biotype', 'gene_source',
                    'note', 'other']
gtf_df = gtf_df[['gene_id', 'gene_biotype', 'note']]
gtf_df['gene_id'] = gtf_df['gene_id'].str.replace("gene_id ", "")
gtf_df['gene_id'] = gtf_df['gene_id'].str.replace('"', '')
gtf_df['gene_biotype'] = gtf_df['gene_biotype'].str.replace("gene_biotype ", "")
gtf_df['gene_biotype'] = gtf_df['gene_biotype'].str.replace('"', '')
gtf_df['note'] = gtf_df['note'].str.replace("note ", "")
gtf_df['note'] = gtf_df['note'].str.replace('"', '')

# Merge gtf gene_biotype and note to tpm_df
tpm_df = tpm_df.merge(gtf_df, how='left', left_on='gene_id', right_on='gene_id')

# Create an average TPM column per sample type
sample_type = ['WT', 'KO', 'MLP1_IP']

for i, sample in enumerate(sample_type):
    triplicates = list(tpm_df.filter(like=sample).columns)
    tpm_df[sample+'_avg_TPM'] = tpm_df[triplicates].mean(axis=1)


tpm_df.to_csv(snakemake.output.tpm_df_biotype, index=False, sep="\t")
