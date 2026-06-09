import polars as pl
import polars_bio as pb
import json


# This file contains STAR-aligned reads from simulated BEERS2 data
# we will check for reads that originated in pseudogenes and count how many are multimappers or nonmappers
R1 = "/home/thobr/nonuniform_impact/data/data/all_bias/beers/results/S1_L1_R1.fastq"
R2 = "/home/thobr/nonuniform_impact/data/data/all_bias/beers/results/S1_L1_R2.fastq"
OUT_R1 = "processed/BEERS_no_pseudogene.R1.fastq"
OUT_R2 = "processed/BEERS_no_pseudogene.R2.fastq"

config = json.load(open("config.json", "rt"))

GRCm38_annot = pb.scan_gtf(
    config["MOUSE_GRCm38_GTF"],
    attr_fields=[
        "gene_id",
        "gene_name",
        "gene_biotype",
        "transcript_id",
        "transcript_biotype",
        "transcript_version",
        "exon_id",
    ],
)
# Note: we have to collect here and not later to avoid a bug where the second filter() gets dropped
GRCm38_transcript_annot = GRCm38_annot.filter(pl.col("type") == "transcript").collect()

pseudogenes = set(
    GRCm38_transcript_annot.filter(
        pl.col("gene_biotype").str.contains("pseudogene")
    ).select("transcript_id")["transcript_id"]
)

# R1
no_pseudogenes = (
    pb.scan_fastq(R1)
    .with_columns(transcript_id=pl.col("name").str.split("_").list.get(1))
    .filter(~pl.col("transcript_id").is_in(pseudogenes))
)
pb.sink_fastq(
    no_pseudogenes,
    OUT_R1,
)

# R2
no_pseudogenes = (
    pb.scan_fastq(R2)
    .with_columns(transcript_id=pl.col("name").str.split("_").list.get(1))
    .filter(~pl.col("transcript_id").is_in(pseudogenes))
)
pb.sink_fastq(
    no_pseudogenes,
    OUT_R2,
)
