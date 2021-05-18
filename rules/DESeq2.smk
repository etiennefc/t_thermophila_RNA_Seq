import os

rule DESeq2_genes:
    """ Differential expression per gene for the different conditions. The
        samples in design.tsv must match exactly the values (not the keys) of
        the 'datasets' dictionary located in the config.json file."""
    input:
        counts = os.path.join(config["path"]["coco_merge"], "merged_counts.tsv"),
        samples = "data/design.tsv"
    output:
        results = directory("results/DESeq2/genes"),
    log:
        "logs/DESeq2/genes.log"
    conda:
        "../envs/deseq2_new.yaml"
    script:
        "../scripts/DESeq2_genes.R"
