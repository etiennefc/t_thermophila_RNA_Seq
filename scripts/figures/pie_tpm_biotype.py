#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv(snakemake.input.tpm_df, sep='\t')
biotype_colors_dict = snakemake.params.biotype_colors
# Generate a pie chart per average sample type (WT, KO, MLP1_IP) of total
# abundance per gene biotype
tpm_list = ft.tpm_list(df, ['WT', 'KO', 'MLP1_IP'], '_avg_TPM', [' tRNA',
            ' rRNA', ' protein_coding'], 'gene_biotype')
tpm_percent = ft.percent_count(tpm_list)

ft.pie_multiple(tpm_percent, ['tRNA', 'rRNA', 'protein_coding'],
                list(biotype_colors_dict.values()),
                ['Average of WT samples', 'Average of MLP1 KO samples',
                'Average of MLP1 IP samples'],
                'Proportion of the total abundance (in TPM) per biotype',
                snakemake.output.pie)
