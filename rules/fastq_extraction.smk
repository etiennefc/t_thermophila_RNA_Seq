import os

rule extract_fastq_reads:
    """ Extract reads from untrimmed fastq corresponding to specific
        isotype_anticodon tRNAs (based on their unique sequence and
        allowing 2 mismatches across the tRNA sequence)."""
    input:
        tRNA_sequences = config['path']['tRNA_unique_sequences'],
        fastq = "data/references/fastq/{sample_ids}.fastq.gz"
    output:
        fasta_dir = directory(os.path.join(config['path']['extract_fastq_reads'], '{sample_ids}/'))
    params:
        agrep_path = "git_repos/agrep/agrep"
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/extract_fastq_reads.py"

rule extract_read_type:
    """ Extract the number of read per type of tRNA (mature ending with CCA or
        premature ending with T{1,10}) for each tRNA isotype_anticodon and per
        condition."""
    input:
        fasta_dir = os.path.join(config['path']['extract_fastq_reads'], '{sample_ids}/'),
        tRNA_sequences = config['path']['tRNA_unique_sequences']
    output:
        read_type_dir = directory(os.path.join(config['path']['extract_read_type'], '{sample_ids}/'))
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/extract_read_type.py"

rule merge_isotypes_per_sample:
    """ Merge together all tRNA isotype_anticodon reads tables generated by
        extract_read_type per sample."""
    input:
        read_type_dir = rules.extract_read_type.output.read_type_dir
    output:
        empty_log = 'logs/{sample_ids}_merge_isotypes_per_sample.log'
    params:
        read_type_df = 'results/extract_read_type/{sample_ids}/merged_isotypes.tsv'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_isotypes_per_sample.py"

rule normalize_read_tables:
    """ For each sample, normalize the number of reads per read_type by the
        total number of reads associated to that isotype_anticodon. Return a
        normalized table per sample."""
    input:
        read_type_df = 'results/extract_read_type/{sample_ids}/merged_isotypes.tsv'
    output:
        empty_log = 'logs/{sample_ids}_normalize_read_tables.log'
    params:
        normalized_df = 'results/extract_read_type/{sample_ids}/merged_isotypes_normalized.tsv'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/normalize_read_tables.py"

rule average_normalized_tables:
    """  Create an average normalized table per condition (WT, KO and MLP1 IP)."""
    input:
        normalized_dfs = expand('results/extract_read_type/{sample_ids}/merged_isotypes_normalized.tsv', sample_ids=config['sample_ids'])
    output:
        WT_df = os.path.join(config['path']['extract_read_type'], 'WT_avg_normalized_df.tsv'),
        KO_df = os.path.join(config['path']['extract_read_type'], 'KO_avg_normalized_df.tsv'),
        MLP1_IP_df = os.path.join(config['path']['extract_read_type'], 'MLP1_IP_avg_normalized_df.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/average_normalized_tables.py"