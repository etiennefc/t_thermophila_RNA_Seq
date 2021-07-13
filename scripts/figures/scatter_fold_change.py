#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

wt_ko_df = pd.read_csv(snakemake.input.deseq_wt_ko, sep='\t')
wt_ip_df = pd.read_csv(snakemake.input.deseq_wt_ip, sep='\t')
biotype_df = pd.read_csv(snakemake.input.biotype_df, sep='\t')
biotype_colors_dict = snakemake.params.biotype_colors

# Merge biotype and avg_TPM columns to the 2 DESeq dfs
wt_ko_df.columns = ['gene_id', 'log2FoldChange_KO', 'pvalue', 'padj', '-log_padj']
wt_ko_df = wt_ko_df[['gene_id', 'log2FoldChange_KO']].merge(biotype_df[['gene_name',
                    'gene_biotype', 'WT_avg_TPM', 'KO_avg_TPM', 'MLP1_IP_avg_TPM']],
                            how='left', left_on='gene_id', right_on='gene_name')
wt_ko_df = wt_ko_df.drop(columns='gene_id')

wt_ip_df.columns = ['gene_id', 'log2FoldChange_IP', 'pvalue', 'padj', '-log_padj']
df = wt_ip_df[['gene_id', 'log2FoldChange_IP']].merge(wt_ko_df,
                            how='left', left_on='gene_id', right_on='gene_name')

# Select RNAs with >1 TPM in the average samples of either WT, KO or MLP1 IP samples
df = df.loc[(df[['WT_avg_TPM', 'KO_avg_TPM', 'MLP1_IP_avg_TPM']] > 1).any(axis=1)]
df.to_csv(snakemake.output.output_df, sep='\t', index=False)

# Create scatter plot of the log2 FC of MLP1 KO (relatively to WT) vs log2 FC of MLP1 IP (relatively to WT)
ft.scatter(df, 'log2FoldChange_IP', 'log2FoldChange_KO', 'gene_biotype',
            'Log2 of the fold change between'+'\n'+'MLP1 IP and WT samples',
            'Log2 of the fold change between'+'\n'+'MLP1 KO and WT samples',
            'Comparison of the log2 fold changes of expressed RNAs'+'\n'+'between MLP1 KO and MLP1 IP (relatively to WT)',
            biotype_colors_dict, snakemake.output.scatter, alpha=0.7)
