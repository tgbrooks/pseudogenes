import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    import json
    import simple_gtf
    import polars_bio as pb
    import lets_plot as lp

    return json, lp, mo, pb, pl


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Sanity checks on the pseudogenes
    """)
    return


@app.cell
def _(pl):
    gene_mapping_raw = pl.read_csv("results/human.parent_gene_mapping.txt", separator="\t")
    gene_mapping_raw
    return (gene_mapping_raw,)


@app.cell
def _(json, pb):
    config = json.load(open("config.json", "rt"))

    annot = pb.read_gtf(config["HUMAN_GTF"], attr_fields= [
        "gene_id",
        "gene_name",
        "gene_biotype",
        "transcript_id",
        "transcript_biotype",
        "transcript_version",
        "exon_id",
    ])
    return annot, config


@app.cell
def _(annot):
    annot
    return


@app.cell
def _(annot):
    transcript_annot = annot.filter(type="transcript")
    gene_annot = annot.filter(type="gene")
    return gene_annot, transcript_annot


@app.cell
def _(gene_annot, gene_mapping_raw, transcript_annot):
    # Attach gene names and biotypes
    gene_mapping = gene_mapping_raw.join(
        transcript_annot.select(
            pseudogene_id = "transcript_id",
            pseudogene_gene_id = "gene_id",
            pseudogene_biotype = "transcript_biotype",
        ),
        how="left",
        on="pseudogene_id",
    ).join(
        gene_annot.select(
            pseudogene_gene_id = "gene_id",
            pseudogene_name = "gene_name",
        ),
        how="left",
        on="pseudogene_gene_id",
    ).join(
        gene_annot.select(
            gene_id = "gene_id",
            gene_name = "gene_name",
            gene_biotype = "gene_biotype",
        ),
        how="left",
        on="gene_id",
    )
    gene_mapping
    return (gene_mapping,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Compare our blastn-inferred parents to what we expect from the gene names

    Most pseudogenes have a gene name of the format "{gene name}P{number}" so we look for these and check if their parent name is "{gene name}" as expected.

    In practice, we see most genes match. Some pseudogenes do not have the expected naming scheme and others names would imply they map to a gene that does not exist in the annotation. We don't consider these as true mismatches since the name was not informative. Many of the mismatch genes are going to a different gene in the same family, see table below. Another source of mismatches are unitary pseudogenes that are annotated in Ensembl twice, once as a unitary pseudogene and once as a lncRNA or occasionally as protein-coding.
    """)
    return


@app.cell
def _(annot, gene_mapping, pl):
    _mapped = gene_mapping.select(
        "pseudogene_gene_id",
        "pseudogene_name",
        "gene_name",
        inferred_gene_name = pl.when(
            pl.col("pseudogene_name").str.contains("^[A-Z0-9]+P[0-9]*$"))
        .then(pl.col("pseudogene_name").str.extract("^([A-Z0-9]+)P[0-9]*$", group_index=1))
        .otherwise(pl.lit(None)),
        pident = "pident",
        pseudogene_biotype = "pseudogene_biotype",
        gene_biotype = "gene_biotype",
        #length = "length",
    ).with_columns(
        inferred_gene_name_exists = pl.col("inferred_gene_name").is_in(set(annot['gene_name'].unique()))
    )
    _class = _mapped.select(
        pl.when(pl.col("inferred_gene_name").is_null())
        .then(pl.lit("1. name not formatted"))
        .when(~pl.col("inferred_gene_name_exists"))
        .then(pl.lit("2. no such gene name"))
        .when(pl.col("inferred_gene_name") == pl.col("gene_name"))
        .then(pl.lit("3. matching gene name (success)"))
        .otherwise(pl.lit("4. nonmatching gene name (failure)"))
        .alias("class")
    )['class'].value_counts().sort("class")
    print(_class)
    _mapped.filter(pl.col("inferred_gene_name") != pl.col("gene_name"))
    return


@app.cell
def _(annot, pl):
    annot.filter(pl.col("feature").is_in(['gene', 'transcript']), gene_id = "ENSG00000279170")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Sequence similarity to parent gene
    """)
    return


@app.cell
def _(gene_mapping, pl):
    gene_mapping.select(
        "pseudogene_id",
        "pseudogene_gene_id", 
        "gene_id",
        pl.col("length"),
        pl.col("slen"),
        pl.col("qlen"),
        length_frac = pl.col("length") / pl.col("slen"),
        pident = "pident",
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Alignment of simulated reads

    We use BEERS2 simulated reads from a mouse dataset to assess the potential impact of pseudogene expression on other gene types, particularly protein coding genes.
    """)
    return


@app.cell
def _(lp, pl):
    alignment_results = pl.read_csv("results/BEERS_STAR_alignment_by_transcript.txt", separator="\t")
    _data = (
        alignment_results
            .filter(
                pl.col("n_unique") + pl.col("n_multi") + pl.col("n_unmapped") > 100
            ).drop_nulls(subset=["gene_id"])
            .with_columns(
                pl.when(pl.col('is_pseudogene')).then(pl.lit("pseudogene")).otherwise(pl.lit("gene")).alias("class")
            )
        .filter(pl.len().over("gene_biotype") >= 5)
    )
    _p2 = lp.ggplot(_data, lp.aes(x="gene_biotype", y = "p_multi", color="is_pseudogene")) + lp.geom_boxplot()
    _p2 + lp.ggsize(width=800,height=550) + lp.theme(axis_text_x = lp.element_text(angle = 45, hjust = 1, vjust = 1)) \
        + lp.labs(y = "percent multimappers", x="")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Here we have checked which simulated reads map to each gene, so we have a mapping from source transcripts to which genes it mapped to.
    """)
    return


@app.cell
def _(config, pb):
    GRCm38_annot = pb.read_gtf(
        config["MOUSE_GRCm38_GTF"],
        attr_fields= [
            "gene_id",
            "gene_name",
            "gene_biotype",
            "transcript_id",
            "transcript_biotype",
            "transcript_version",
            "exon_id",
        ],
    )
    GRCm38_transcript_annot = GRCm38_annot.filter(type ="transcript")
    GRCm38_gene_annot = GRCm38_annot.filter(type ="gene")
    return GRCm38_gene_annot, GRCm38_transcript_annot


@app.cell
def _(GRCm38_gene_annot, GRCm38_transcript_annot, pl):
    mapping = (
        pl.read_csv("processed/BEERS2_transcripts_to_STAR_mapped_genes.txt", separator="\t")
        .join(
            GRCm38_transcript_annot.select(
                source_transcript = "transcript_id",
                source_gene_id = "gene_id",
            ),
            on = "source_transcript",
        )
        .join(
            GRCm38_gene_annot.select(
                source_gene_id = "gene_id",
                source_gene_name = "gene_name",
                source_gene_biotype = "gene_biotype",
            ),
            on="source_gene_id"
        )
        .join(
            GRCm38_gene_annot.select(
                mapped_gene_id = "gene_id",
                mapped_gene_name = "gene_name",
                mapped_gene_biotype = "gene_biotype",
            ),
            on="mapped_gene_id"
        )
    )
    from_pseudogenes = mapping.filter(pl.col("source_gene_biotype").str.contains("pseudogene"))
    to_pseudogenes = mapping.filter(pl.col("mapped_gene_biotype").str.contains("pseudogene"))
    #from_pseudogenes.filter(mapped_gene_biotype = "protein_coding").sort("count", descending=True)
    return from_pseudogenes, mapping, to_pseudogenes


@app.cell
def _(from_pseudogenes, mapping, mo, pl):
    total_reads = 65958848 # actually just mapped reads from samtools view -c -F 260 -q 255
    number_of_misplaced_reads = from_pseudogenes.filter(mapped_gene_biotype = "protein_coding", multimapper=False).sort("count", descending=True).select(pl.col("count").sum()).item()
    number_of_misplaced_reads
    print(f"Fraction of reads misplaced: {number_of_misplaced_reads / total_reads:0.2%}")

    misplaced_by_gene = mapping.filter(multimapper=False, mapped_gene_biotype="protein_coding").group_by(
        "mapped_gene_id", "mapped_gene_name",
    ).agg(
        total_reads = pl.col("count").sum(),
        misplaced_pseudogene_reads = pl.col("count").filter(pl.col("source_gene_biotype").str.contains("pseudogene")).sum()
    ).with_columns(
        percent_misplaced = pl.col("misplaced_pseudogene_reads") / pl.col("total_reads")
    )
    mo.ui.table(
        misplaced_by_gene.sort("percent_misplaced", descending=True).filter(pl.col('total_reads') > 100),
        label = "Misplaced reads mapping to protein coding genes, arising from pseudogenes",
        format_mapping = {"percent_misplaced": "{:0.2%}".format},
        selection = None,
    )
    return (total_reads,)


@app.cell
def _(mapping, mo, pl, to_pseudogenes, total_reads):
    number_of_misplaced_reads2 = to_pseudogenes.filter(source_gene_biotype = "protein_coding", multimapper=False).sort("count", descending=True).select(pl.col("count").sum()).item()
    print(f"Fraction of reads misplaced: {number_of_misplaced_reads2 / total_reads:0.2%}")

    misplaced_by_gene2 = mapping.filter(
        pl.col("mapped_gene_biotype").str.contains("pseudogene"),
        multimapper=False, 
    ).group_by(
        "mapped_gene_id", "mapped_gene_name",
    ).agg(
        total_reads = pl.col("count").sum(),
        misplaced_protein_coding_reads = pl.col("count").filter(source_gene_biotype = "protein_coding").sum()
    ).with_columns(
        percent_misplaced = pl.col("misplaced_protein_coding_reads") / pl.col("total_reads")
    )
    mo.ui.table(
        misplaced_by_gene2.sort("percent_misplaced", descending=True).filter(pl.col('total_reads') > 100),
        label = "Misplaced reads mapping to pseudogenes, arising from protein coding genes",
        format_mapping = {"percent_misplaced": "{:0.2%}".format},
        selection = None,
    )
    return


@app.cell
def _(mapping):
    #mapping.filter(mapped_gene_name = "Gapdh", multimapper=False)
    mapping.filter(source_gene_name = "Gm7336")
    return


@app.cell
def _(GRCm38_gene_annot):
    GRCm38_gene_annot.filter(gene_name = "Gm7336")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Salmon quants from BEERS simulated data
    Here we compare Salmon quants using BEERS data but we have removed all reads originating in pseudogenes.
    """)
    return


@app.cell
def _(GRCm38_gene_annot, pl):
    salmon_quant_schema =  {"Name": pl.Utf8, "Length": pl.Float64, "EffectiveLength": pl.Float64, "TPM": pl.Float64, "NumReads": pl.Float64}
    gene_quants = pl.concat([
        pl.read_csv("processed/BEERS_no_pseudogene/salmon/quant.genes.sf", separator="\t", schema=salmon_quant_schema)
            .with_columns(type=pl.lit("no_pseudogenes")),
        pl.read_csv("processed/BEERS_all_transcripts/salmon/quant.genes.sf", separator="\t", schema=salmon_quant_schema)
            .with_columns(type=pl.lit("all_transcripts")),
    ]).rename({"Name": "gene_id"}).join(
        GRCm38_gene_annot,
        on="gene_id",
    )
    return (gene_quants,)


@app.cell
def _(gene_quants):
    gene_quants
    return


app._unparsable_cell(
    r"""
    _wide = gene_quants.pivot(
        on = "type",
        index=["gene_id", "gene_biotype"],
        values="TPM",
    )
    _p1 = lp.ggplot(_wide, lp.aes(y="no_pseudogenes", x="all_transcripts")) + lp.geom_pointdensity() + lp.scale_x_log10() + lp.scale_y_log10() + lp.scale_color_viridis() + lp.ggtitle("TPM before/after removing pseudogene reads") + lp.coord_fixed()
    _p2 = lp.ggplot(_wide.filter(~pl.col("gene_biotype").str.contains("pseudogene")), lp.aes(y="no_pseudogenes", x="all_transcripts")) + lp.geom_pointdensity() + lp.scale_x_log10() + lp.scale_y_log10() + lp.scale_color_viridis() + lp.ggtitle("pseudogenes removed") + lp.coord_fixed()
    lp.gggrid([_p1, _p2]) + lp.ggsize(1000, 500) \
    """,
    name="_"
)


@app.cell
def _(pl):
    true_TPM = pl.read_parquet("/home/thobr/nonuniform_impact/data/beers.true_TPM.parquet").rename({"GeneID": "gene_id", "TranscriptID": "transcript_id", "TPM": "true_TPM"}) \
        .filter(sample = 1) # We only use one sample in this
    true_gene_TPM = true_TPM.group_by("gene_id").agg(pl.col("true_TPM").sum())
    return (true_gene_TPM,)


@app.cell
def _(gene_quants, pl, true_gene_TPM):
    def normalize(c):
        return pl.col(c) / pl.col(c).sum() * 1e6
    
    TPM_wide = gene_quants.pivot(
        on = "type",
        index=["gene_id", "gene_biotype", "chrom"],
        values="TPM",
    ).join(true_gene_TPM, on="gene_id")\
    .filter(
        # these are very high expressed and dominate TPM normalization of other genes
        # We drop MT genes and then renormalize everything
        pl.col("chrom") != "MT",
    ).with_columns(
        true_TPM = normalize("true_TPM"),
        no_pseudogenes = normalize("no_pseudogenes"),
        all_transcripts = normalize("all_transcripts"),
    )
    return (TPM_wide,)


@app.cell
def _(TPM_wide, lp, pl):
    _p1 = lp.ggplot(TPM_wide.filter(~pl.col("gene_biotype").str.contains("pseudogene")), lp.aes(y="all_transcripts", x="true_TPM")) + lp.geom_pointdensity() + lp.scale_x_log10() + lp.scale_y_log10() + lp.scale_color_viridis() + lp.ggtitle("all transcripts") + lp.coord_fixed()
    _p2 = lp.ggplot(TPM_wide.filter(~pl.col("gene_biotype").str.contains("pseudogene")), lp.aes(y="no_pseudogenes", x="true_TPM")) + lp.geom_pointdensity() + lp.scale_x_log10() + lp.scale_y_log10() + lp.scale_color_viridis() + lp.ggtitle("pseudogenes removed") + lp.coord_fixed()
    lp.gggrid([_p1, _p2]) + lp.ggsize(1000, 500)
    return


@app.cell
def _(TPM_wide, mo, pl):
    # Summarize TPM performance with/without pseudogenes
    def MARD(A,B):
        A,B = pl.col(A), pl.col(B)
        return ((A-B).abs() / (A+B)).median()
    mo.vstack([
        "MARD (median absolute relative deviance) comparing true TPM to Salmon TPMs on non-pseudogenes. FASTQ either with all BEERS simulated reads or with just the non-pseudogene reads. NOTE: lower is better.",
        TPM_wide
            .filter(
                ~pl.col("gene_biotype").str.contains("pseudogene"),
                pl.col("true_TPM") > 10, # drop very low-expressed
            )
            .select(
            MARD_all_transcripts = MARD("true_TPM", "all_transcripts"),
            MARD_no_pseudogenes = MARD("true_TPM", "no_pseudogenes"),
        )
    ])
    return


if __name__ == "__main__":
    app.run()
