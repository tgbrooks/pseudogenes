# Pseudogenes investigation

Here, we investigate pseudogenes with an interest on their impact on quantification of (short-read) RNA-seq.
Pseudogenes often have sequence similarity to genes and so are sources of quantifications mistakes which could impact both pseudogene and gene expression values.
Our goal is to establish:

1. The extent to which pseudogenes are expressed.
2. The extent to which gene expression is misclassified as pseudogene expression.
3. Whether modifications of the annotation (removal of problematic pseudogenes) improves quantifications of genes.

Pseudogenes are classified into the following categories:

1. Processed pseudogene: arising from reverse transcription of RNA and therefore not containing introns.
2. Unprocessed pseudogene: arising from gene duplication and so containing intronic sequence.
3. Unitary pseudogene: arising from LoF mutation in a gene which is present species-wide. Has no paralog gene (but has ortholog genes in other species).

Ensembl has a large number of biotypes that relate to pseudogenes, indicating additional properties:

IG_pseudogene
IG_C_pseudogene
IG_J_pseudogene
IG_V_pseudogene
TR_V_pseudogene
TR_J_pseudogene
Mt_tRNA_pseudogene
tRNA_pseudogene
snoRNA_pseudogene
snRNA_pseudogene
scRNA_pseudogene
rRNA_pseudogene
misc_RNA_pseudogene
miRNA_pseudogene
pseudogene
processed_pseudogene
polymorphic_pseudogene
retrotransposed # This one isn't used in the latest annotations, so all pseudogene biotypes have 'pseudogene' in them
transcribed_processed_pseudogene
transcribed_unprocessed_pseudogene
transcribed_unitary_pseudogene
unitary_pseudogene
unprocessed_pseudogene

# TODO:

1. Determine mapping of pseudogenes to their parent paralog gene (if any) by aligning exons to other genes
2. Quantify overlap between pseudogenes and genes: how many variants separate them? How many intronic variants separate their genomic region?
3. Check upstream/downstream of pseudogene / gene for similarity
4. Assess impact on alignment


## Commands to recreate:
``` bash
module load /project/itmatlab/sharedmodules/salmon-v1.9.0
module load samtools/1.20

auv run python scripts/map_to_parent_genes.py --species human
auv run python scripts/check_beers_alignments.py
auv run python scripts/transcript_to_mapped_gene.py
auv run python scripts/strip_pseudogenes_from_bam.py
auv run python util/add_transcript_versions.py "/project/itmatlab/genomes/mouse/GRCm38/Ensembl.v102/Mus_musculus.GRCm38.102.gtf" "processed/Salmon.GRCm38.gene_transcript_map.txt"
auv run python

# Map with Salmon BEERS fastqs with pseudogenes removed or not
export SALMON_GRCm38_INDEX="/project/itmatlab/index/SALMON-1.9.0_indexes/GRCm38.ensemblv102/salmon_index"
export SALMON_ARGS="-l A --softclip --softclipOverhangs -p 6 --gcBias --posBias --seqBias"
bsub -M 36000 -R "rusage [mem=36000] span[hosts=1]" -n 6 -eo "logs/Salmon.no_pseudogenes.err" -oo "logs/Salmon.no_pseudogenes.out" \
    salmon quant -i $SALMON_GRCm38_INDEX -g "processed/Salmon.GRCm38.gene_transcript_map.txt" $SALMON_ARGS \
        -1 "processed/BEERS_no_pseudogene.R1.fastq" -2 "processed/BEERS_no_pseudogene.R2.fastq" \
        -o "processed/BEERS_no_pseudogene/salmon/"

bsub -M 36000 -R "rusage [mem=36000] span[hosts=1]" -n 6 -eo "logs/Salmon.all_transcripts.err" -oo "logs/Salmon.all_transcripts.out" \
    salmon quant -i $SALMON_GRCm38_INDEX -g "processed/Salmon.GRCm38.gene_transcript_map.txt" $SALMON_ARGS \
        -1 "/home/thobr/nonuniform_impact/data/data/all_bias/beers/results/S1_L1_R1.fastq" \
        -2 "/home/thobr/nonuniform_impact/data/data/all_bias/beers/results/S1_L1_R2.fastq" \
        -o "processed/BEERS_all_transcripts/salmon/"

```

# Notebooks

Marimo notebooks have final analyses, figures.

```
auv run marimo edit notebooks/
```
