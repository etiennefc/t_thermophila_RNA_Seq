#!/usr/bin/python3
from snakemake.io import expand
import os
def get_figures_path(config):
    """Return list of figures (and their path) to generate from config"""
    files = []
    #tpm = "{attributes}_tpm_N24_vs_N0.svg"
    #files.extend(expand(os.path.join(config['figure']['scatter'], 'read_diff_IP_vs_WT.svg'), **config))
    #files.extend(expand(os.path.join(config['figure']['scatter'], 'read_diff_KO_vs_WT.svg'), **config))
    files.extend(expand(os.path.join(config['figure']['scatter'], 'read_diff_IP_vs_WT_wo_trna.svg'), **config))
    files.extend(expand(os.path.join(config['figure']['scatter'], 'read_diff_KO_vs_WT_wo_trna.svg'), **config))        
    #files.append(os.path.join(config['figure']['volcano'], 'KO_vs_WT.svg'))
    #files.append(os.path.join(config['figure']['volcano'], 'MLP1_IP_vs_WT.svg'))
    #files.append(os.path.join(config['figure']['volcano'], 'KO_vs_WT_tpm.svg'))
    #files.append(os.path.join(config['figure']['volcano'], 'MLP1_IP_vs_WT_tpm.svg'))
    #files.append(os.path.join(config['figure']['pie'], 'average_TPM_per_sample.svg'))
    #files.append(os.path.join(config['figure']['scatter'], 'log2FC_KO_vs_IP.svg'))


    print(files)
    return files
