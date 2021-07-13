#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv(snakemake.input.df, sep='\t')

# Create scatter plot of the number of reads difference between premature and
# mature tRNAs by comparing MLP1 IP avg sample to WT avg samples
ft.scatter(df, 'MLP1_IP_read_diff', 'WT_read_diff', None,
            'Difference of number of reads between the premature and'+'\n'+'mature form of a tRNA in average MLP1 IP samples',
            'Difference of number of reads between the premature and'+'\n'+'mature form of a tRNA in average WT samples',
            'Comparison of premature vs mature tRNAs read accumulations'+'\n'+'between MLP1 IP and WT average samples)',
            'grey', snakemake.output.ip_vs_wt)

# Create scatter plot of the number of reads difference between premature and
# mature tRNAs by comparing MLP1 KO avg sample to WT avg samples
ft.scatter(df, 'KO_read_diff', 'WT_read_diff', None,
            'Difference of number of reads between the premature and'+'\n'+'mature form of a tRNA in average KO samples',
            'Difference of number of reads between the premature and'+'\n'+'mature form of a tRNA in average WT samples',
            'Comparison of premature vs mature tRNAs read accumulations'+'\n'+'between KO and WT average samples)',
            'grey', snakemake.output.ko_vs_wt)
