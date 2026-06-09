# NOTE: this snakefile doesn't really work! We have issues with python and snakemake
# I needed a newer version of GLIB_C for some dependencies (polars-bio) but that requires
# running it through an apptainer container, which doesn't allow submitting of jobs!
# Trying to run the python stuff in their own container for some reason fought over the .venv
# with the outside python (for snakemake). My resolution is to do it by hand

configfile: "config.json"

rule all:
    input:
        "results/human.parent_gene_mapping.txt",
        #"results/biotype_summary.txt",
        "processed/BEERS_no_pseudogene/salmon/quant.sf",
        "processed/BEERS_all_transcripts/salmon/quant.sf",

rule map_to_parent_genes:
    input:
        gtf = config['HUMAN_GTF'],
        fa = config['HUMAN_FA'],
    output:
        "results/human.parent_gene_mapping.txt",
    resources:
        mem_mb = 36_000,
        threads=6,
    container:
        config['PYTHON_CONTAINER']
    script:
        "uv run python scripts/map_to_parent_genes.py --species human"

rule check_beers_alignments:
    input:
        gtf = config['MOUSE_GRCm38_GTF'],
        bam = "/home/thobr/nonuniform_impact/processed/unbiased/samples/S1/bam/Aligned.sortedByCoord.out.bam",
    output:
        "results/BEERS_STAR_alignment_by_transcript.txt",
    resources:
        mem_mb = 36_000,
        threads = 6,
    container:
        config['PYTHON_CONTAINER']
    script:
        "uv run python scripts/check_beers_alignments.py"

rule transcript_to_mapped_genes:
    input:
        gtf = config['MOUSE_GRCm38_GTF'],
        bam = "/home/thobr/nonuniform_impact/processed/unbiased/samples/S1/bam/Aligned.sortedByCoord.out.bam",
    output:
        "processed/BEERS2_transcripts_to_STAR_mapped_genes.txt"
    resources:
        mem_mb = 36_000,
        threads = 6,
    container:
        config['PYTHON_CONTAINER']
    script:
        "uv run python scripts/transcript_to_mapped_gene.py"

rule strip_pseudogenes_from_fastq:
    input:
        gtf = config['MOUSE_GRCm38_GTF'],
        R1 = "/home/thobr/nonuniform_impact/data/data/all_bias/beers/results/S1_L1_R1.fastq",
        R2 = "/home/thobr/nonuniform_impact/data/data/all_bias/beers/results/S1_L1_R2.fastq",
    output:
        "processed/BEERS_no_pseudogene.R1.fastq"
        "processed/BEERS_no_pseudogene.R2.fastq"
    resources:
        mem_mb = 6_000,
    container:
        config['PYTHON_CONTAINER']
    script:
        "uv run python scripts/strip_pseudogenes_from_bam.py"

rule salmon_transcript_map:
    input:
        gtf = config['MOUSE_GRCm38_GTF']
    output:
        txt = "processed/Salmon.GRCm38.gene_transcript_map.txt"
    localrule: True # so we can use the python version here
    shell:
        "auv run python util/add_transcript_versions.py {input} {output}"

rule salmon_quant_no_pseudogenes:
    input:
        index = config['MOUSE_GRCm38_SALMON_INDEX'],
        gene_map = "processed/Salmon.GRCm38.gene_transcript_map.txt",
        R1 = "processed/BEERS_no_pseudogene.R1.fastq",
        R2 = "processed/BEERS_no_pseudogene.R2.fastq",
    output:
        "processed/BEERS_no_pseudogene/salmon/quant.sf"
    resources:
        mem_mb = 25_000,
        threads = 6,
    params:
        args = "-l A --softclip --softclipOverhangs -p 6 --gcBias --posBias --seqBias",
        out_folder = "processed/BEERS_no_pseudogene/salmon/"
    shell:
        "salmon quant -i {input.index} -g {input.gene_map} {params.args} -1 {input.R1} -2 {input.R2} -o {params.out_folder}"

rule salmon_quant_all_transcripts:
    input:
        index = config['MOUSE_GRCm38_SALMON_INDEX'],
        gene_map = "processed/Salmon.GRCm38.gene_transcript_map.txt",
        R1 = "/home/thobr/nonuniform_impact/data/data/all_bias/beers/results/S1_L1_R1.fastq",
        R2 = "/home/thobr/nonuniform_impact/data/data/all_bias/beers/results/S1_L1_R2.fastq",
    output:
        "processed/BEERS_all_transcripts/salmon/quant.sf"
    resources:
        mem_mb = 25_000,
        threads = 6,
    params:
        args = "-l A --softclip --softclipOverhangs -p 6 --gcBias --posBias --seqBias",
        out_folder = "processed/BEERS_all_transcripts/salmon/"
    shell:
        "salmon quant -i {input.index} -g {input.gene_map} {params.args} -1 {input.R1} -2 {input.R2} -o {params.out_folder}"
