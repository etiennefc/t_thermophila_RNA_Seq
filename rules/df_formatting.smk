import os

include: "DESeq2.smk"

rule add_gene_biotype:
    """ Add gene biotype to all genes quantified by CoCo cc and an average
        abundance column per sample type."""
    input:
        gtf = config['path']['complete_gtf'],
        tpm_df = os.path.join(config['path']['coco_merge'], "merged_tpm.tsv")
    output:
        tpm_df_biotype = config['path']['tpm_df_biotype']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get_gene_biotype.py"

rule format_deseq_output:
    """ Format DESeq2 output to keep only genes with non-null values for log2
        fold change and p-value."""
    input:
        WT_vs_KO = os.path.join(rules.DESeq2_genes.output.results, 'partial_knockout-wild_type.csv'),
        WT_vs_IP = os.path.join(rules.DESeq2_genes.output.results, 'MLP1_IP-wild_type.csv')
    output:
        WT_vs_KO_clean = os.path.join(config['path']['deseq_output'],
                        'partial_knockout-wild_type_v2.csv'),
        WT_vs_IP_clean = os.path.join(config['path']['deseq_output'],
                        'MLP1_IP-wild_type_v2.csv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/format_deseq_output.py"
