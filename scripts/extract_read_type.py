#!/usr/bin/python3
import pandas as pd
import subprocess as sp

tRNA = pd.read_csv(snakemake.input.tRNA_sequences)
input_fasta_dir = snakemake.input.fasta_dir
output_dir = snakemake.output.read_type_dir

# Create dict of isotype_anticodon as keys and their fishing sequence(s) in a list as values
d = {}
for isotype in list(pd.unique(tRNA['Isotype'])):
    seqs = list(tRNA[tRNA['Isotype'] == isotype].sequence.values)
    d[isotype] = seqs

sp.call(f'mkdir -p {output_dir}', shell=True)
for k, v in d.items():
    sp.call(f'touch {output_dir}/{k}.tsv', shell=True)
    # Mature tRNA reads ending with CCA or CCAN (where N can be any nt)
    sp.call(f'grep -E "CCA.?$" {input_fasta_dir}/{k}.fa | wc -l >> {output_dir}/{k}.tsv', shell=True)
    # Premature tRNAs ending with 1 to 10 T
    for i in range(10):
        j = i + 1
        sp.call('grep -E "[^T]T{'+str(j)+'}'+f'$" {input_fasta_dir}/{k}.fa | wc -l >> {output_dir}/{k}.tsv', shell=True)
    # Other tRNA reads not fitting in previous categories
    sp.call('grep -vE "(CCA.?$)|(T{1,10}$)|(^>)" '+f'{input_fasta_dir}/{k}.fa | wc -l >> {output_dir}/{k}.tsv', shell=True)
