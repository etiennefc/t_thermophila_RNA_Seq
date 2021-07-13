#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

wt_ko_df = pd.read_csv(snakemake.input.deseq_wt_ko, sep='\t')
wt_ip_df = pd.read_csv(snakemake.input.deseq_wt_ip, sep='\t')
tpm_biotype_df = pd.read_csv(snakemake.input.tpm_biotype_df, sep='\t')
biotype_colors_dict = snakemake.params.biotype_colors
tpm_biotype_colors_dict = snakemake.params.biotype_tpm_colors

# Merge biotype and avg_TPM columns to the 2 DESeq dfs
wt_ko_df = wt_ko_df.merge(tpm_biotype_df[['gene_name', 'gene_biotype', 'WT_avg_TPM', 'KO_avg_TPM']],
                            how='left', left_on='gene_id', right_on='gene_name')

wt_ip_df = wt_ip_df.merge(tpm_biotype_df[['gene_name', 'gene_biotype', 'WT_avg_TPM', 'MLP1_IP_avg_TPM']],
                            how='left', left_on='gene_id', right_on='gene_name')


# Generate a volcano plot for the two conditions compared to WT (hue: gene_biotype)
ft.volcano(wt_ko_df, 'log2FoldChange', '-log_padj', 'gene_biotype',
            'Log2 of the fold change between'+'\n'+'MLP1 KO and WT samples',
            '-log10(adjusted p-value)',
            'Abundance of all RNAs'+'\n'+'in MLP1 KO compared to WT',
            biotype_colors_dict, snakemake.output.volcano_wt_ko, alpha=0.7)

ft.volcano(wt_ip_df, 'log2FoldChange', '-log_padj', 'gene_biotype',
            'Log2 of the fold change between'+'\n'+'MLP1 IP and WT samples',
            '-log10(adjusted p-value)',
            'Abundance of all RNAs'+'\n'+'in MLP1 IP compared to WT',
            biotype_colors_dict, snakemake.output.volcano_wt_ip, alpha=0.7)

#Create a tpm_biotype hue column (based on if the RNA is > 1TPM in at least of the two average sample that are compared)
#wt_ko_df.loc[(wt_ko_df['gene_biotype'] == ' rRNA') & (wt_ko_df[['WT_avg_TPM', 'KO_avg_TPM']] > 1).any(axis=1), 'tpm_biotype_ko'] = 'rRNA_expressed'
#wt_ko_df.loc[(wt_ko_df['gene_biotype'] == ' rRNA') & (wt_ko_df[['WT_avg_TPM', 'KO_avg_TPM']] <= 1).all(axis=1), 'tpm_biotype_ko'] = 'rRNA_not_expressed'
#wt_ko_df.loc[(wt_ko_df['gene_biotype'] == ' tRNA') & (wt_ko_df[['WT_avg_TPM', 'KO_avg_TPM']] > 1).any(axis=1), 'tpm_biotype_ko'] = 'tRNA_expressed'
#wt_ko_df.loc[(wt_ko_df['gene_biotype'] == ' tRNA') & (wt_ko_df[['WT_avg_TPM', 'KO_avg_TPM']] <= 1).all(axis=1), 'tpm_biotype_ko'] = 'tRNA_not_expressed'
#wt_ko_df.loc[(wt_ko_df['gene_biotype'] == ' protein_coding') & (wt_ko_df[['WT_avg_TPM', 'KO_avg_TPM']] > 1).any(axis=1), 'tpm_biotype_ko'] = 'pc_expressed'
#wt_ko_df.loc[(wt_ko_df['gene_biotype'] == ' protein_coding') & (wt_ko_df[['WT_avg_TPM', 'KO_avg_TPM']] <= 1).all(axis=1), 'tpm_biotype_ko'] = 'pc_not_expressed'

#wt_ip_df.loc[(wt_ip_df['gene_biotype'] == ' rRNA') & (wt_ip_df[['WT_avg_TPM', 'MLP1_IP_avg_TPM']] > 1).any(axis=1), 'tpm_biotype_ip'] = 'rRNA_expressed'
#wt_ip_df.loc[(wt_ip_df['gene_biotype'] == ' rRNA') & (wt_ip_df[['WT_avg_TPM', 'MLP1_IP_avg_TPM']] <= 1).all(axis=1), 'tpm_biotype_ip'] = 'rRNA_not_expressed'
#wt_ip_df.loc[(wt_ip_df['gene_biotype'] == ' tRNA') & (wt_ip_df[['WT_avg_TPM', 'MLP1_IP_avg_TPM']] > 1).any(axis=1), 'tpm_biotype_ip'] = 'tRNA_expressed'
#wt_ip_df.loc[(wt_ip_df['gene_biotype'] == ' tRNA') & (wt_ip_df[['WT_avg_TPM', 'MLP1_IP_avg_TPM']] <= 1).all(axis=1), 'tpm_biotype_ip'] = 'tRNA_not_expressed'
#wt_ip_df.loc[(wt_ip_df['gene_biotype'] == ' protein_coding') & (wt_ip_df[['WT_avg_TPM', 'MLP1_IP_avg_TPM']] > 1).any(axis=1), 'tpm_biotype_ip'] = 'pc_expressed'
#wt_ip_df.loc[(wt_ip_df['gene_biotype'] == ' protein_coding') & (wt_ip_df[['WT_avg_TPM', 'MLP1_IP_avg_TPM']] <= 1).all(axis=1), 'tpm_biotype_ip'] = 'pc_not_expressed'


# Generate a volcano plot for the two conditions compared to WT only for expressed RNAs (hue: gene_biotype)
wt_ko_df = wt_ko_df.loc[(wt_ko_df[['WT_avg_TPM', 'KO_avg_TPM']] > 1).any(axis=1)]
ft.volcano(wt_ko_df, 'log2FoldChange', '-log_padj', 'gene_biotype',
            'Log2 of the fold change between'+'\n'+'MLP1 KO and WT samples',
            '-log10(adjusted p-value)',
            'Abundance of expressed RNAs'+'\n'+'in MLP1 KO compared to WT',
            biotype_colors_dict, snakemake.output.volcano_wt_ko_tpm, alpha=0.7)
wt_ip_df = wt_ip_df.loc[(wt_ip_df[['WT_avg_TPM', 'MLP1_IP_avg_TPM']] > 1).any(axis=1)]
ft.volcano(wt_ip_df, 'log2FoldChange', '-log_padj', 'gene_biotype',
            'Log2 of the fold change between'+'\n'+'MLP1 IP and WT samples',
            '-log10(adjusted p-value)',
            'Abundance of expressed RNAs'+'\n'+'in MLP1 IP compared to WT',
            biotype_colors_dict, snakemake.output.volcano_wt_ip_tpm, alpha=0.7)
