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
retrotransposed
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
