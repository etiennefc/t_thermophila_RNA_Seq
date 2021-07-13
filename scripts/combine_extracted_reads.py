#!/usr/bin/python3
import pandas as pd


""" Generate a dataframe where each column corresponds to the number of reads
    associated to either premature tRNA or mature tRNA in all the samples (18
    columns). Each line corresponds to a tRNA."""
path = snakemake.params.path
sample_files = ['KO_1', 'KO_2', 'KO_3', 'MLP1_IP_1', 'MLP1_IP_2',
                'MLP1_IP_3', 'WT_1', 'WT_2', 'WT_3']
output = snakemake.output.merged_reads_df
# For each sample, create a dataframe of 2 columns (1 for the number of mature
# tRNA reads, one for the number of premature tRNA reads)
dfs = []
for name in sample_files:
    with open(path+name+'.txt', 'r') as file:
        lines = file.read().split('\n')[0:-1]

        # Select tRNA ids (each 5 lines starting from the first line) into a list
        tRNA_id = [id.replace('>', '') for id in lines[0::5]]

        # Select number of mature tRNA reads (each 5 lines starting from the 3rd line) into a list
        mature = [int(read) for read in lines[2::5]]

        # Select number of premature tRNA reads (each 5 lines starting from the 4th line) into a list
        premature = [int(read) for read in lines[4::5]]

        # Create dictionary
        d = {name+'_mature': mature, name+'_premature': premature}

        # Create dataframe
        df = pd.DataFrame(d, columns=[name+'_mature', name+'_premature'], index=tRNA_id)
        dfs.append(df)

# Concat horizontally all df into one df and create a gene_id column by resetting the index
final_df = pd.concat(dfs, axis=1)
final_df = final_df.reset_index().rename(columns={'index':'gene_id'})
final_df.to_csv(output, index=False, sep='\t')
