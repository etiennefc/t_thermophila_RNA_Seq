#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv(snakemake.input.df, sep='\t')

# Generate a scatter plot of the abundance of expressed intronic snoRNAs at N0 and N24
# with different hues (snoRNA type, correlation of abundance with host gene and host gene biotype)
ft.scatter(df, 'N0_avg', 'N24_avg', snakemake.wildcards.attributes,
            'Average abundance at N0 (TPM)', 'Average abundance at N24 (TPM)', '',
            snakemake.params.hue, snakemake.output.scatter, hue_order=snakemake.params.hue.keys())


# Generate a volcano plot of the log2(fold change between N24 vs N0) and -log10(adjusted pval)
# with different hues (snoRNA type, correlation of abundance with host gene and host gene biotype)
ft.volcano(df, 'log2FoldChange', 'log_padj', snakemake.wildcards.attributes,
            'Log2 of the fold change between N24 and N0', '-log10(adjusted p-value)', '',
            snakemake.params.hue, snakemake.output.volcano, hue_order=snakemake.params.hue.keys())
