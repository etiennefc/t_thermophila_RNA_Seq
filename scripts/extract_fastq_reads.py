#!/usr/bin/python3
import pandas as pd
import subprocess as sp
import regex, gzip, time

tRNA = pd.read_csv(snakemake.input.tRNA_sequences)
input_fastq = snakemake.input.fastq
output_dir = snakemake.output.fasta_dir
agrep_path = snakemake.params.agrep_path
sample_id = snakemake.wildcards.sample_ids

# Create dict of isotype_anticodon as keys and their fishing sequence(s) in a list as values
d = {}
for isotype in list(pd.unique(tRNA['Isotype'])):
    seqs = list(tRNA[tRNA['Isotype'] == isotype].sequence.values)
    d[isotype] = seqs

# For each isotype_anticodon, create a file containing all reads containing its possible corresponding fishing sequence(s) with <=2 mismatch allowed
sp.call(f'mkdir -p {output_dir}', shell=True)
temp_fastq = f'temp_{sample_id}.fastq'
sp.call(f'gunzip -c {input_fastq} > {temp_fastq}', shell=True)
for k, v in d.items():
    if len(v) == 1:  # if only one fishing sequence per isotype_anticodon
        seq = v[0]
        sp.call(f'echo {k} {seq}', shell=True)
        sp.call(f'{agrep_path} -2 {seq} {temp_fastq} | head -n -1 > {output_dir}/{k}.fa', shell=True)
        sp.call("""gawk -i inplace 'BEGIN{j=0} {j+=1; print ">read_"j"\\n"$0}' """+f"""{output_dir}/{k}.fa""", shell=True)  # add read id for all reads
    else:  # if multiple fishing sequence per isotype_anticodon
        len_v = len(v) - 1
        for i, seq in enumerate(v):
            sp.call(f'echo {k} {seq}', shell=True)
            sp.call(f'{agrep_path} -2 {seq} {temp_fastq} | head -n -1 >> {output_dir}/{k}.fa', shell=True)
            sp.call(f'wc -l {output_dir}/{k}.fa', shell=True)
            if i == len_v:  # for last iteration, add read id for all reads that match all sequences in v
                sp.call("""gawk -i inplace 'BEGIN{j=0} {j+=1; print ">read_"j"\\n"$0}' """+f"""{output_dir}/{k}.fa""", shell=True)

sp.call(f'rm {temp_fastq}', shell=True)




# Python fuzzy regex search implementation (WAY longer than with bash agrep only!!)
"""
# For each isotype_anticodon, create a file containing all reads containing its possible corresponding fishing sequence(s) with <=2 mismatch allowed
sp.call(f'mkdir -p {output_dir}', shell=True)
for k, v in d.items():
    sp.call(f'touch {output_dir}/{k}.fa', shell=True)
    if len(v) == 1:  # if only one fishing sequence per isotype_anticodon
        seq = v[0]
        sp.call(f'echo {k} {seq}', shell=True)
        start = time.time()
        with gzip.open(input_fastq, 'rb') as f:  # read fastq and search for reads containing that sequence
            i = 1
            for line in f:
                line = line.decode('utf-8')
                m = regex.findall("("+seq+"){s<=2}", line)  # allow max two substitution across the sequence
                if len(m) >= 1:  # if sequence is present in the read
                    l = line.strip('\n')
                    sp.call(f'echo ">"read_{k}_{i} >> {output_dir}/{k}.fa', shell=True)  # add fasta identifier for this read
                    sp.call(f'echo {l} >> {output_dir}/{k}.fa', shell=True) # add matched read
                    i += 1
        end = time.time()
        print(end - start)
    else:  # if multiple fishing sequence per isotype_anticodon
        adapted_seq = ['('+j+')' for j in v]  # Add '()' around each sequences
        regex_or = '{s<=2}|'.join(adapted_seq)+'{s<=2}' # Create regex of multiple "or" ex: "(ATTTC){s<=2}|(CGCGCT){s<=2}|(TTTAGG){s<=2}"
        sp.call(f'echo {k}', shell=True)
        print(len(v))
        start = time.time()
        with gzip.open(input_fastq, 'rb') as f:  # read fastq and search for reads containing that sequence
            i = 1
            for line in f:
                line = line.decode('utf-8')
                m = regex.findall(regex_or, line)  # allow max two substitution across each possible sequence
                if len(m) >= 1:  # if sequence is present in the read
                    l = line.strip('\n')
                    sp.call(f'echo ">"read_{k}_{i} >> {output_dir}/{k}.fa', shell=True)
                    sp.call(f'echo {l} >> {output_dir}/{k}.fa', shell=True)
                    i += 1
        end = time.time()
        print(end - start)
"""
