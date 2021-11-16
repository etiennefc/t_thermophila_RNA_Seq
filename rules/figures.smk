import os

include: "df_formatting.smk"
include: "bam_extraction.smk"
include: "fastq_extraction.smk"

rule pie_chart:
    """ Create a pie chart for all average sample of the cumulative abundance
        (in TPM) per biotype. """
    input:
        tpm_df = rules.add_gene_biotype.output.tpm_df_biotype
    output:
        pie = os.path.join(config['figure']['pie'], 'average_TPM_per_sample.svg')
    params:
        biotype_colors = config['colors']['gene_biotype']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/figures/pie_tpm_biotype.py"

rule volcano:
    """ Create a volcano plot of the two comparisons (MLP1_IP vs WT and KO vs
        WT (fold change relative to WT))."""
    input:
        deseq_wt_ko = rules.format_deseq_output.output.WT_vs_KO_clean,
        deseq_wt_ip = rules.format_deseq_output.output.WT_vs_IP_clean,
        tpm_biotype_df = rules.add_gene_biotype.output.tpm_df_biotype
    output:
        volcano_wt_ko = os.path.join(config['figure']['volcano'], 'KO_vs_WT.svg'),
        volcano_wt_ip = os.path.join(config['figure']['volcano'], 'MLP1_IP_vs_WT.svg'),
        volcano_wt_ko_tpm = os.path.join(config['figure']['volcano'], 'KO_vs_WT_tpm.svg'),
        volcano_wt_ip_tpm = os.path.join(config['figure']['volcano'], 'MLP1_IP_vs_WT_tpm.svg')
    params:
        biotype_colors = config['colors']['gene_biotype'],
        biotype_tpm_colors = config['colors']['tpm_biotype']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/figures/volcano.py"

rule scatter_fold_change:
    """ Create a scatter plot of the log2 fold change of either MLP1 KO or IP vs
        WT samples."""
    input:
        deseq_wt_ko = rules.format_deseq_output.output.WT_vs_KO_clean,
        deseq_wt_ip = rules.format_deseq_output.output.WT_vs_IP_clean,
        biotype_df = rules.add_gene_biotype.output.tpm_df_biotype
    output:
        scatter = os.path.join(config['figure']['scatter'], 'log2FC_KO_vs_IP.svg'),
        output_df = config['path']['merged_deseq_fold_change']
    params:
        biotype_colors = config['colors']['gene_biotype']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/figures/scatter_fold_change.py"

rule scatter_read_diff:
    """ Create a scatter plot of the difference of number of premature vs mature
        tRNA reads by comparing either IP vs WT or KO vs WT read differences."""
    input:
        df = rules.mature_premature_diff.output.avg_diff_df
    output:
        #ip_vs_wt = os.path.join(config['figure']['scatter'], 'read_diff_IP_vs_WT.svg'),
        #ko_vs_wt = os.path.join(config['figure']['scatter'], 'read_diff_KO_vs_WT.svg')
        ip_vs_wt = os.path.join(config['figure']['scatter'], 'read_diff_IP_vs_WT_wo_trna.svg'),
        ko_vs_wt = os.path.join(config['figure']['scatter'], 'read_diff_KO_vs_WT_wo_trna.svg')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/figures/scatter_read_diff.py"

rule heatmap_read_type:
    """ Create one heatmap per average condition to show the number of reads per
        read type and tRNA isotype."""
    input:
        WT_df = rules.average_normalized_tables.output.WT_df,
        KO_df = rules.average_normalized_tables.output.KO_df,
        MLP1_IP_df = rules.average_normalized_tables.output.MLP1_IP_df
    output:
        heatmap_wt = os.path.join(config['figure']['heatmap'], 'WT_read_type_per_isotype.svg'),
        heatmap_ko = os.path.join(config['figure']['heatmap'], 'KO_read_type_per_isotype.svg'),
        heatmap_ip = os.path.join(config['figure']['heatmap'], 'MLP1_IP_read_type_per_isotype.svg')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/figures/heatmap_read_type.py"

rule heatmap_read_type_ratio_wt_ip:
    """ Create one heatmap to show the ratio of number of reads per
        read type and tRNA isotype between MLP1 IP and WT average samples."""
    input:
        WT_df = rules.average_normalized_tables.output.WT_df,
        MLP1_IP_df = rules.average_normalized_tables.output.MLP1_IP_df
    output:
        heatmap = os.path.join(config['figure']['heatmap'], 'read_type_ratio_wt_ip.svg')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/figures/heatmap_read_type_ratio_wt_ip.py"
