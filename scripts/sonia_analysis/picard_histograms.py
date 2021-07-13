#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import sys

"""
Script producing histograms from columns 'insert_size' and 'All_Reads.fr_count' contained in 
Picard_insert_size_metrics.txt files. Each of these file must be in a separate folder named with the sample name. You 
must change the sample_lst according to your sample names.
"""

sample_lst = ['KO_1', 'KO_2', 'KO_3', 'WT_1', 'WT_2', 'WT_3', 'MLP1_IP_1', 'MLP1_IP_2', 'MLP1_IP_3']

def txt_to_df(text_path):
    """Extract useful data from txt file and return dataframe."""
    with open(str(text_path+'/Picard_insert_size_metrics.txt'), 'r+') as file:
        while not 'All_Reads.fr_count' in next(file):
            pass
        temp = []
        for line in file:
            data = line.strip().split()
            temp.append(data)

        df = pd.DataFrame(temp, columns=['insert_size', 'all_read_counts'])
        df = df.dropna()
        print(df)
        return df

def histogram(df, sample_name, figure_output_path):
    """Generate histogram from given dataframe."""
    x = [int(item) for item in list(df['insert_size'])]
    y = [int(item) for item in list(df['all_read_counts'])]
    df = pd.DataFrame({'insert_size':x, 'all_read_counts':y})
    ax = df.plot.bar(x='insert_size', y='all_read_counts', color='blue')
    ax.set_xlabel('Insert size', fontsize=12)
    ax.set_ylabel('Total number of inserts', fontsize=12)
    plt.title(sample_name, fontsize=15)
    
    ax.set_xticklabels(x, fontsize=10, rotation=0)
    legend = ax.legend()
    legend.remove()

    every_nth = 50
    for n, label in enumerate(ax.xaxis.get_ticklabels()):
        if n % every_nth != 0:
            label.set_visible(False)

    plt.savefig(figure_output_path, dpi=600, bbox_inches='tight')


def main(sample_common_path, sample_list=sample_lst):
    """Iterate through sample list to generate a histogram for each sample (stored in the same folder than the txt
        file)"""
    for i, sample in enumerate(sample_list):
        df = txt_to_df(str(sample_common_path)+sample)
        histogram(df, sample, str(sample_common_path)+'/'+sample+'/'+sample+'.png')





main('/home/fcouture/scratch/t_thermophila_rna_seq/results/Picard/')


























































