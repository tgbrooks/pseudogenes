import polars as pl
import pysam
import simple_gtf
import json


# This file contains STAR-aligned reads from simulated BEERS2 data
# we will check for reads that originated in pseudogenes and count how many are multimappers or nonmappers
BAM_FILE = "/home/thobr/nonuniform_impact/processed/unbiased/samples/S1/bam/Aligned.sortedByCoord.out.bam"

##################################################
# Read annotations
config = json.load(open("config.json", "rt"))
annot_raw = simple_gtf.read_gtf(config["MOUSE_GTF"])

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


##################################################
# Go through the reads in the bam file
bam = pysam.AlignmentFile(BAM_FILE, "rb")
by_transcript = {}
for i, read in enumerate(bam):
    # We only consider one alignment per fragment
    if read.is_secondary:
        continue
    if not read.is_read1:
        continue

    transcript_id = read.qname.split(":")[6].split("_")[1]
    if transcript_id not in by_transcript:
        by_transcript[transcript_id] = {"n_unmapped": 0, "n_multi": 0, "n_unique": 0}
    res = by_transcript[transcript_id]
    if not read.is_mapped:
        res["n_unmapped"] += 1
    else:
        num_loci_mapped = read.get_tag("NH")
        if num_loci_mapped > 1:
            res["n_multi"] += 1
        else:
            res["n_unique"] += 1
    if i % 1_000_000 == 0:
        print(".", end="", flush=True)
print("done")

results = pl.DataFrame(
    [{"transcript_id": k, **res} for k, res in by_transcript.items()]
).join(
    annot.filter(feature="transcript").select(
        "gene_id", "gene_biotype", "transcript_id"
    ),
    how="left",
    on="transcript_id",
)

gene_results = (
    results.group_by("gene_id", "gene_biotype")
    .agg(
        pl.col("n_unmapped").sum(),
        pl.col("n_multi").sum(),
        pl.col("n_unique").sum(),
    )
    .with_columns(
        is_pseudogene=pl.col("gene_biotype").cast(str).str.contains("pseudogene"),
        p_unmapped=pl.col("n_unmapped")
        / (pl.col("n_unmapped") + pl.col("n_multi") + pl.col("n_unique")),
        p_multi=pl.col("n_multi")
        / (pl.col("n_unmapped") + pl.col("n_multi") + pl.col("n_unique")),
        p_unique=pl.col("n_unique")
        / (pl.col("n_unmapped") + pl.col("n_multi") + pl.col("n_unique")),
    )
)

gene_results.write_csv("results/BEERS_STAR_alignment_by_transcript.txt", separator="\t")
