import os

rule get_tRNA_sequences:
    """ Extract tRNA sequences from gtf file (coordinates) and genome fasta file."""
    input:
        gtf = config['path']['complete_gtf'],
        genome_fasta = config['path']['reference_genome']
    output:
        tRNA_sequences = config['path']['tRNA_sequences']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get_tRNA_sequences.py"

rule position_to_gene_id:
    """ Modify fasta file of tRNA sequences (from get_tRNA_sequences) so that
        the position lines (starting with '>' plus chr_00...) are changed to add
        the tRNA id after the position (with a '<' between the position and the
        gene_id)."""
    input:
        original_fasta = rules.get_tRNA_sequences.output.tRNA_sequences,
        gtf = config['path']['complete_gtf']
    output:
        modified_fasta = config['path']['tRNA_sequences_id']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/position_to_gene_id.py"

rule extract_reads:
    """ For each tRNA in each sample, extract the number of reads associated to
        either mature tRNAs (i.e. ending with CCA) or premature tRNAs (i.e.
        ending with a track of 0-5 T). For mature tRNAs, we take the last 50 nt
        (or 49 nt since sometimes this nt is removed before adding the CCA) and
        look for reads containing that sequence plus a CCA. For premature tRNAs,
        we take the last 50 nt and look for reads containing that sequence plus
        a serie of 1-5 T.
        ***Does not work with WT samples and KO_2... Obliged to run them locally one by one..."""
    input:
        tRNA_sequences = rules.position_to_gene_id.output.modified_fasta,
        bam = "results/star/{id}/Aligned.sortedByCoord.out.bam"
    output:
        number_reads = os.path.join(config['path']['extract_reads'], "{id}.txt")
    params:
        script = "scripts/extract_reads.sh"
    conda:
        "../envs/coco.yaml"
    shell:
        "bash {params.script} {input.tRNA_sequences} {input.bam} {output.number_reads}"

rule merge_reads:
    """ Merge dataframe of mature and premature tRNA reads for each sample
        (output from extract_reads) into one dataframe."""
    params:
        #path = os.path.join(config['path']['extract_reads'], 'Original/'),
        path = config['path']['extract_reads_wo_trna_seq']
    output:
        merged_reads_df = config['path']['merged_reads_df']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/combine_extracted_reads.py"

rule mature_premature_diff:
    """ Get the difference between the number of premature tRNA reads and mature
        tRNA reads for all average samples."""
    input:
        df = rules.merge_reads.output.merged_reads_df
    output:
        avg_diff_df = config['path']['mature_premature_diff']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/mature_premature_diff.py"

rule extract_reads_wo_trna_seq:
    """ For each tRNA in each sample, extract the number of reads associated to either mature
        tRNAs (ending with CCA) or premature tRNAs (ending with a track of 1-5 T) We don't look 
        for the sequence of the tRNA before the CCA or serie of T."""
    input:
        df = rules.position_to_gene_id.output.modified_fasta,
        bam = "results/star/{id}/Aligned.sortedByCoord.out.bam" 
    output:
        extracted_reads = os.path.join(config['path']['extract_reads_wo_trna_seq'], "{id}.txt")
    shell:
        "./scripts/extract_reads_wo_trna_seq.sh {input.df} {input.bam} {output.extracted_reads}"

rule merge_reads_wo_trna_seq:
    """ Merge dataframe of mature and premature tRNA reads for each sample
        (output from extract_reads_wo_trna_seq) into one dataframe."""
    params:
        path = config['path']['extract_reads_wo_trna_seq']
    output:
        merged_reads_df = config['path']['merged_reads_wo_trna_df']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/combine_extracted_reads.py"

rule mature_premature_diff_wo_trna_seq:
    """ Get the difference between the number of premature tRNA reads and mature
        tRNA reads for all average samples."""
    input:
        df = rules.merge_reads_wo_trna_seq.output.merged_reads_df
    output:
        avg_diff_df = config['path']['mature_premature_diff_wo_trna']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/mature_premature_diff.py"

