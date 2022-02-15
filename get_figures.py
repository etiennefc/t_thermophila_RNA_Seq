#!/usr/bin/python3
from snakemake.io import expand
import os
def get_figures_path(config):
    """Return list of figures (and their path) to generate from config"""
    files = []
    files.append(os.path.join(config['figure']['heatmap'], 'WT_read_type_per_isotype.svg'))
    files.append(os.path.join(config['figure']['heatmap'], 'KO_read_type_per_isotype.svg'))
    files.append(os.path.join(config['figure']['heatmap'], 'MLP1_IP_read_type_per_isotype.svg'))
    files.append(os.path.join(config['figure']['heatmap'], 'read_type_ratio_wt_ip.svg'))

    print(files)
    return files
