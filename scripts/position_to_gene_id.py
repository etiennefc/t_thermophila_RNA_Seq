#!/usr/bin/python3
import pandas as pd

""" Add gene_id after position in tRNA sequence file."""

original_fasta = snakemake.input.original_fasta
cols = ['chr', 'source', 'feature', 'start', 'end', 'score1', 'strand', 'score2', 'features']
gtf = pd.read_csv(snakemake.input.gtf, sep='\t', names=cols)
tRNA_gtf = gtf[gtf['feature'] == 'gene']
tRNA_gtf = tRNA_gtf[tRNA_gtf['features'].str.contains('gene_biotype "tRNA"')]  # select tRNA genes only
tRNA_gtf['gene_id'] = tRNA_gtf['features'].str[:14].str.replace('gene_id "', '')  # create gene_id column

# Create position column
tRNA_gtf['position'] = '>' + tRNA_gtf['chr'] + ':' + tRNA_gtf['start'].astype(str) + '-' + tRNA_gtf['end'].astype(str) + '\n'

# Create dictionary of positions (ex: >chr_096:88458-88530) as a key and gene_id as the value (ex: 80630)
ref_dict = tRNA_gtf.set_index('position').to_dict()['gene_id']

# Add gene_id after position in tRNA sequence fasta file
with open(original_fasta, 'r') as file:
    with open(snakemake.output.modified_fasta, 'w') as output:
        for line in file:
            if line.startswith('>'):
                id = ref_dict[line]
                position = line[1:].replace('\n', '')
                output.write('>'+position+'<'+id+'\n')
            else:
                output.write(line)
