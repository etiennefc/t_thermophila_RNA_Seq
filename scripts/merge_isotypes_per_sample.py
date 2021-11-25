#!/usr/bin/python3
import os
import pandas as pd

input_dir = snakemake.input.read_type_dir
output_df = snakemake.output.read_type_df
col_list = ['mature_CCA', 'premature_1T', 'premature_2T', 'premature_3T',
            'premature_4T', 'premature_5T', 'premature_6T', 'premature_7T',
            'premature_8T', 'premature_9T', 'premature_10T', 'other_reads']

# Add column names to df (see col_list)
dfs = []
for file in os.listdir(input_dir):
    isotype_anticodon = file.split('.')[0]
    if '.tsv' in file:
        df = pd.read_csv(f'{input_dir}/{file}', names=[isotype_anticodon])
        df['cols'] = col_list # temp column that will be transposed in multiple columns (one new col per row)
        df = df.set_index('cols')
        df = df.T
        dfs.append(df)

# Concat vertically all isotype dfs into one df (one per sample)
final_df = pd.concat(dfs)
final_df = final_df.reset_index()
final_df = final_df.rename(columns = {'index': 'isotype_anticodon'})
final_df = final_df.sort_values(by=['isotype_anticodon'])
final_df.to_csv(output_df, index=False, sep='\t')


