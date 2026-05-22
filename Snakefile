configfile: "config.json"

rule all:
    "results/parent_gene_mapping.txt",
    "results/biotype_summary.txt",

rule map_to_parent_genes:
    input:
        gtf = HUMAN_GTF,
        fa = HUMAN_FA,
    output:
        "results/parent_gene_mapping.txt",
        "results/biotype_summary.txt",
    resources:
        mem_mb = 36_000,
        threads=6,
    script:
        "scripts/map_to_parent_genes.py"
