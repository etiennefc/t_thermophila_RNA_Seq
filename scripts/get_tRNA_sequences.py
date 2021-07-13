#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp

""" Extract sequence of all tRNAs from genome fasta file and coordinates from
    the gtf file. We excluded 4 tRNAs since they are out of the range of the
    fasta file of the genome."""

cols = ['chr', 'source', 'feature', 'start', 'end', 'score1', 'strand', 'score2', 'features']
gtf = pd.read_csv(snakemake.input.gtf, sep='\t', names=cols)
fasta = snakemake.input.genome_fasta
out_of_range_tRNAs = ['80234', '80406', '80407', '80540']

# Create bed file of all tRNAs in the annotation
tRNA_gtf = gtf[gtf['feature'] == 'gene']
tRNA_gtf = tRNA_gtf[tRNA_gtf['features'].str.contains('gene_biotype "tRNA"')]
tRNA_gtf['gene_id'] = tRNA_gtf['features'].str[:14].str.replace('gene_id "', '')  # select gene id from features column
tRNA_bed = tRNA_gtf[['chr', 'start', 'end', 'gene_id', 'score1', 'strand', 'source', 'feature', 'score2', 'features']]
tRNA_bed = tRNA_bed.loc[~tRNA_bed['gene_id'].isin(out_of_range_tRNAs)]  # exclude out of range tRNAs
tRNA_bed.to_csv('temp.bed', index=False, sep='\t', header=False)
tRNA_bed = BedTool('temp.bed')
seq = tRNA_bed.sequence(fi=fasta)

# The seqfn attribute points to the fasta of sequences of seq
with open(seq.seqfn, 'r') as file:
    with open(snakemake.output.tRNA_sequences, 'w') as output:
        for line in file:
            output.write(line)


sp.call('rm temp.bed', shell=True)
