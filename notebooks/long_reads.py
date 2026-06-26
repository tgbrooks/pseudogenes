import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Pseudogenes in long read data
    We use a dataset we generated with Shaon Sengupta on mouse PerDKO mice versus WT that were sequenced with both short and long reads.
    Moreover, the short reads were sequenced with both PolyA and a rRNA removal method.
    """)
    return


@app.cell
def _():
    import marimo as mo
    import polars as pl
    import polars_bio as pb
    import lets_plot as lp

    return mo, pb, pl


@app.cell
def _(pb):
    # Load the annotation
    GTF_FILE = "/project/itmatlab/genomes/mouse/GRCm39/Ensembl.v114/Mus_musculus.GRCm39.114.gtf.gz"
    annot = pb.read_gtf(GTF_FILE, attr_fields=[
        "gene_id",
        "gene_name",
        "gene_biotype",
        "transcript_id",
        "transcript_biotype",
        "transcript_version",
        "exon_id",
        "tag",
    ])
    return (annot,)


@app.cell
def _(annot, pl):
    # map pseudogenes to parent genes
    gene_mapping_raw = pl.read_csv("results/mouse.parent_gene_mapping.txt", separator="\t")
    # Attach gene names and biotypes
    gene_mapping = gene_mapping_raw.join(
        annot.filter(type="transcript").select(
            pseudogene_id = "transcript_id",
            pseudogene_gene_id = "gene_id",
            pseudogene_biotype = "transcript_biotype",
        ),
        how="left",
        on="pseudogene_id",
    ).join(
        annot.filter(type="gene").select(
            pseudogene_gene_id = "gene_id",
            pseudogene_name = "gene_name",
        ),
        how="left",
        on="pseudogene_gene_id",
    ).join(
        annot.filter(type="gene").select(
            gene_id = "gene_id",
            gene_name = "gene_name",
            gene_biotype = "gene_biotype",
        ),
        how="left",
        on="gene_id",
    )
    gene_mapping
    return (gene_mapping,)


@app.cell
def _(annot):
    annot.filter(type="gene", gene_name="Gapdh")
    return


@app.cell
def _(pl):
    # Load the bam
    sample_id = "6-182"
    SAMPLE_INFO_FILE = "/project/itmatlab/FIRST_LONG_READ_EXPERIMENT/METADATA/samples.csv"
    sample_info = pl.read_csv(SAMPLE_INFO_FILE).with_columns(sample_id="sample")
    long_read_id = dict(zip(sample_info['sample_id'].str.replace("\\.","-"), sample_info['enumeration']))
    long_read_bam_file = f"/project/itmatlab/FIRST_LONG_READ_EXPERIMENT/RAW_DATA_PROCESSING/ALIGNMENT/flnc-{long_read_id[sample_id]}.bam"
    long_read_bam_file
    return (long_read_bam_file,)


@app.cell
def _(pl):
    long_read_counts_raw = pl.read_csv("processed/short_vs_long_reads/quants.6-128.long_reads.txt", separator="\t", comment_prefix="#")
    long_read_counts_raw = long_read_counts_raw.rename({long_read_counts_raw.columns[-1]: "counts"})
    long_read_counts_raw
    return (long_read_counts_raw,)


@app.cell
def _(pl):
    short_read_counts_raw = pl.read_csv("processed/short_vs_long_reads/quants.6-128.short_reads.txt", separator="\t", comment_prefix="#")
    short_read_counts_raw = short_read_counts_raw.rename({short_read_counts_raw.columns[-1]: "counts"})
    short_read_counts_raw
    return (short_read_counts_raw,)


@app.cell
def _(annot, long_read_counts_raw, short_read_counts_raw):
    counts = long_read_counts_raw.select(gene_id = "Geneid", long_read_counts = "counts", gene_length="Length") \
        .join(
            short_read_counts_raw.select(gene_id = "Geneid", short_read_counts = "counts"),
            on = "gene_id",
        ).join(
            annot.filter(type="gene"),
            on = "gene_id",
        )
    counts
    return (counts,)


@app.cell
def _(counts, gene_mapping, pl):
    # Get counts by protein-coding gene and by its psuedogenes
    _df = counts.join(gene_mapping.select("pseudogene_gene_id", parent_gene_id = "gene_id"), left_on="gene_id", right_on="pseudogene_gene_id", how="left")
    gene_pseudogene_counts = counts.filter(gene_biotype="protein_coding")\
        .select('gene_id', 'gene_name', 'long_read_counts', 'short_read_counts') \
        .join(
            _df.group_by("parent_gene_id")
                .agg(
                    pseudogene_long_read_counts = pl.col("long_read_counts").sum(),
                    pseudogene_short_read_counts = pl.col("short_read_counts").sum(),
                    lr_main_pseudogene = pl.col("gene_id").filter(pl.col("long_read_counts") == pl.col("long_read_counts").max()),
                    sr_main_pseudogene = pl.col("gene_id").filter(pl.col("short_read_counts") == pl.col("short_read_counts").max()),
                ),
            left_on = "gene_id",
            right_on="parent_gene_id",
            how="left",
        ).fill_null(0)
    gene_pseudogene_counts
    return


@app.cell
def _(annot):
    annot.filter(type="gene", gene_id="ENSMUSG00000078965")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Compare sequences at some pseudogenes vs gene
    Now that we see long and short reads differing substantially in some genes (e.g., Gapdh), we want to know if the long reads really are able to distinguish between these.
    Alternatively, it could come down to somewhat arbitrary aligner decisions: STAR prefers annotated junctions over non-spliced reads and therefore could choose a gene over a processed pseudogene without any base-level evidence. Likewise, minimap2 (pbmm2 is what was used) does the opposite by default (in fact it looks like it wasn't provided annotations at all).
    """)
    return


@app.cell
def _(counts, gene_mapping, mo):
    _df = counts.join(
        gene_mapping.select("pseudogene_gene_id", parent_gene_id = "gene_id", parent_gene_name="gene_name"),
        left_on="gene_id",
        right_on="pseudogene_gene_id",
        how="left"
    )
    pseudogenes =  {
        f"{name}: {gene_id} (counts: {counts})": gene_id
        for
        gene_id, parent_gene_id, name, counts in
        _df.sort("long_read_counts", descending=True)
            .select("gene_id", "parent_gene_id", "parent_gene_name", "long_read_counts")
            .drop_nulls()
            .iter_rows()
    }
    gene_selector = mo.ui.dropdown(
        options = pseudogenes,
        searchable=True,
        value=[k for k,v in  pseudogenes.items() if v == "ENSMUSG00000078965"][0]
    )
    return (gene_selector,)


@app.cell
def _(annot, gene_mapping, gene_selector, pl):
    selected_gene_id = gene_selector.value
    _is_primary = pl.col("tag").str.contains("gencode_primary")
    selected_transcript_id = annot.filter(_is_primary, gene_id = selected_gene_id, type="transcript")['transcript_id'][0]
    selected_parent_gene_id = gene_mapping.filter(pseudogene_gene_id = selected_gene_id)['gene_id'][0]
    selected_parent_transcript_id = annot.filter(_is_primary, gene_id = selected_parent_gene_id, type="transcript")['transcript_id'][0]
    return selected_parent_transcript_id, selected_transcript_id


@app.cell
def _(pb, pl):
    TRANSCRIPTOME_FASTAS = ["/project/itmatlab/genomes/mouse/GRCm39/Ensembl.v114/Mus_musculus.GRCm39.cdna.all.fa.gz", "/project/itmatlab/genomes/mouse/GRCm39/Ensembl.v114/Mus_musculus.GRCm39.ncrna.fa.gz"]
    transcriptome = pl.concat([pb.read_fasta(fasta) for fasta in TRANSCRIPTOME_FASTAS]).with_columns(transcript_id=pl.col("name").str.split(".").list.get(0))
    return (transcriptome,)


@app.cell
def _(selected_parent_transcript_id, selected_transcript_id, transcriptome):
    selected_seq = transcriptome.filter(transcript_id = selected_transcript_id)['sequence'][0]
    selected_parent_seq = transcriptome.filter(transcript_id = selected_parent_transcript_id)['sequence'][0]
    print(f"PSEUDOGENE: {selected_seq}")
    print(f"      GENE: {selected_parent_seq}")
    return (selected_parent_seq,)


@app.cell
def _(annot, selected_transcript_id):
    # Here we read one read from that region in the bam file
    # And align it to the gene sequence
    def reverse_complement(seq):
        return seq[::-1].translate(str.maketrans({"A":"T", "C":"G", "G": "C", "T":"A"}))
    def cond_rev_complement(seq, must_reverse):
        if must_reverse:
            return reverse_complement(seq)
        else:
            return seq
    pseudogene_chrom, pseudogene_start, pseudogene_end, pseudogene_strand = next(annot.filter(transcript_id = selected_transcript_id).select("chrom", "start", "end", "strand").iter_rows())
    return (
        pseudogene_chrom,
        pseudogene_end,
        pseudogene_start,
        pseudogene_strand,
        reverse_complement,
    )


@app.cell
def _():
    GRCm39_FASTA = "/project/itmatlab/genomes/mouse/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
    import gzip
    def _():
        genome = {}
        seqs = []
        curr_chrom = None
        for line in gzip.open(GRCm39_FASTA, "rt"):
            if line.startswith(">"):
                if curr_chrom:
                    genome[curr_chrom] = ''.join(seqs)
                    seqs = []
                curr_chrom = line.removeprefix(">").split(" ")[0]
            else:
                seqs.append(line.strip())
        genome[curr_chrom] = ''.join(seqs)
        return genome
    genome = _()
    return (genome,)


@app.cell
def _(
    genome,
    pseudogene_chrom,
    pseudogene_end,
    pseudogene_start,
    pseudogene_strand,
    reverse_complement,
):
    PAD = 2000
    pseudogene_genome_seq = genome[pseudogene_chrom][pseudogene_start-PAD:pseudogene_end+PAD]
    if pseudogene_strand == "-":
       pseudogene_genome_seq = reverse_complement(pseudogene_genome_seq) 
    return (pseudogene_genome_seq,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Align the gene and pseudogene transcripts to each other
    """)
    return


@app.cell
def _():
    import Bio.Align

    # Initialize the aligner
    aligner = Bio.Align.PairwiseAligner()
    aligner.mismatch_score = -1
    aligner.gap_score = -2
    aligner.mode = 'local'
    return (aligner,)


@app.cell
def _(mo):
    read_selector = mo.ui.slider(start=1, stop=100, label="read number: ")
    return (read_selector,)


@app.cell
def _(
    long_read_bam_file,
    pb,
    pl,
    pseudogene_chrom,
    pseudogene_end,
    pseudogene_start,
):
    nearby_reads = (
        pb.scan_bam(long_read_bam_file)
        .filter(
            pl.col("chrom") == pseudogene_chrom,
            pl.col("start") <= pseudogene_end,
            pl.col("end") >= pseudogene_start,
        )
        .select("sequence", "flags", "cigar")
        .head(1000) # Don't take too many, just in case
        .collect()
    )
    return (nearby_reads,)


@app.cell
def _(pseudogene_strand):
    pseudogene_strand
    #nearby_reads.select("flags")
    return


@app.cell
def _(nearby_reads, pl, pseudogene_strand, read_selector, reverse_complement):
    # Find a read that overlaps the pseudogene
    example_read, flags, cigar = next(
            nearby_reads
            .filter(
                # same strand as pseudogene
                (pl.col('flags') & 16 != 0) == (pseudogene_strand == "-"),
                #  no introns - often these just jump over our pseudogene anyway
                #~pl.col("cigar").str.contains("N"),
            )
            .head(read_selector.value).tail(1)
            .iter_rows()
    )
    read_strand = "-" if flags & 16 else "+" # not really needed
    if pseudogene_strand == "-":
        # read sequence in BAM is always on the + strand
        # so we move it to where the pseudogene is.
        # Doesn't actually depend upon the read strand.
        example_read = reverse_complement(example_read)

    return (example_read,)


@app.cell
def _(
    example_read,
    long_read_bam_file,
    pseudogene_chrom,
    pseudogene_end,
    pseudogene_start,
    read_selector,
):
    len(example_read)
    print(f'samtools view {long_read_bam_file} {pseudogene_chrom}:{pseudogene_start}-{pseudogene_end} | head -n {read_selector.value} | tail -n 1')
    print(example_read)
    return


@app.cell(hide_code=True)
def _(gene_selector, mo, read_selector):
    mo.vstack([gene_selector, read_selector])
    return


@app.cell
def _(aligner, example_read, mo, pseudogene_genome_seq, selected_parent_seq):
    # Alignment of pseudogene GENOME sequence to the parent transcript
    from alignment_viewer_widget import AlignmentViewer
    _a1 = aligner.align(selected_parent_seq, pseudogene_genome_seq)[0]
    _a2 = aligner.align(pseudogene_genome_seq, example_read)[0]
    mo.vstack([
        mo.ui.anywidget(AlignmentViewer([_a1, _a2], names=["parent", "pseudogene genome", "read"]))
    ])
    return


if __name__ == "__main__":
    app.run()
