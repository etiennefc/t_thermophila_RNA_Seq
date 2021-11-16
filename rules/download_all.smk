import os


rule download_genome:
    """Download the reference genome (fasta file) used for this analysis."""
    output:
        genome = config['path']['reference_genome']
    params:
        link = config['download']['genome']
    shell:
        "wget --quiet -O {output.genome} {params.link}"


rule download_annotation:
    """Download the annotations (gff files) used for this analysis."""
    output:
        genome_gff = config['path']['genome_gff'],
        tRNA_gff = config['path']['tRNA_gff'],
        rRNA_5s_gff = config['path']['rRNA_5s_gff']
    params:
        link_genome_gff = config['download']['genome_gff'],
        link_tRNA_gff = config['download']['tRNA_gff'],
        link_5s_gff = config['download']['rRNA_5s_gff']
    shell:
        "wget --quiet -O {output.genome_gff} {params.link_genome_gff} && "
        "wget --quiet -O {output.tRNA_gff} {params.link_tRNA_gff} && "
        "wget --quiet -O {output.rRNA_5s_gff} {params.link_5s_gff}"

rule download_coco_git:
    """Download git repository of CoCo."""
    output:
        git_coco_folder = directory('git_repos/coco')
    params:
        git_coco_link = config['path']['coco_git_link']
    conda:
        '../envs/git.yaml'
    shell:
        'mkdir -p {output.git_coco_folder} '
        '&& git clone {params.git_coco_link} {output.git_coco_folder}'

rule download_agrep:
    """ Download git repository of agrep, a bash tool to to approximate (fuzzy)
        grep searches (ex: allow N mismatches in sequence search)."""
    output:
        git_agrep_folder = directory('git_repos/agrep')
    params:
        git_agrep_link = config['path']['agrep_git_link']
    conda:
        '../envs/git.yaml'
    shell:
        'mkdir -p {output.git_agrep_folder} && '
        'git clone {params.git_agrep_link} {output.git_agrep_folder} && '
        'cd git_repos/agrep && make'    
