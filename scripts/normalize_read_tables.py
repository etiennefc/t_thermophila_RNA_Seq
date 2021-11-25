#!/usr/bin/python3
import pandas as pd

df = pd.read_csv(snakemake.input.read_type_df, sep='\t')
output_df = snakemake.output.normalized_df
col_list = ['mature_CCA', 'premature_1T', 'premature_2T', 'premature_3T',
            'premature_4T', 'premature_5T', 'premature_6T', 'premature_7T',
            'premature_8T', 'premature_9T', 'premature_10T', 'other_reads']

# Normalize by total number of reads per isotype (i.e via the total of the row)
df['sum'] = df.sum(axis=1)
for i, col in enumerate(col_list):
    df[col+'_norm'] = df[col] / df['sum']

df.to_csv(output_df, index=False, sep='\t')





