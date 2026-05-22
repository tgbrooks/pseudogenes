import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    import json
    import simple_gtf

    return json, mo, pl, simple_gtf


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Sanity checks on the pseudogenes
    """)
    return


@app.cell
def _(pl):
    gene_mapping_raw = pl.read_csv("results/parent_gene_mapping.txt", separator="\t")
    gene_mapping_raw
    return (gene_mapping_raw,)


@app.cell
def _(json, pl, simple_gtf):
    config = json.load(open("config.json", "rt"))

    annot_raw = simple_gtf.read_gtf(config["HUMAN_GTF"])

    # Annotations are read into lists but many types only ever have one entry, so we collapse them
    singular_columns = [
        "gene_id",
        "gene_name",
        "gene_biotype",
        "transcript_id",
        "transcript_biotype",
        "transcript_version",
        "exon_id",
    ]
    for col in singular_columns:
        assert annot_raw.select(m=pl.col(col).list.len().max())["m"][0] == 1
    annot = annot_raw.with_columns(**{c: pl.col(c).explode() for c in singular_columns})
    return (annot,)


@app.cell
def _(annot):
    annot
    return


@app.cell
def _(annot):
    transcript_annot = annot.filter(feature="transcript")
    gene_annot = annot.filter(feature="gene")
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
    annot.filter(pl.col("feature").is_in(['gene', 'transcript']), gene_name = "CD99P1")
    return


if __name__ == "__main__":
    app.run()
