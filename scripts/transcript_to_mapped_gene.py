import polars as pl
import polars_bio as pb
from util.cigar import cigar_to_ref_blocks
import json


# This file contains STAR-aligned reads from simulated BEERS2 data
# we will check for reads that originated in pseudogenes and count how many are multimappers or nonmappers
BAM_FILE = "/home/thobr/nonuniform_impact/processed/unbiased/samples/S1/bam/Aligned.sortedByCoord.out.bam"

##################################################
# Read annotations
config = json.load(open("config.json", "rt"))
annot = pb.read_gtf(
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

exons = (
    annot.lazy()
    .filter(type="exon")
    .join(
        annot.lazy().filter(type="transcript").select("transcript_id", "gene_id"),
        how="left",
        on="transcript_id",
    )
    .select("chrom", "start", "end", pl.col("gene_id").cast(pl.Categorical))
    .sort("chrom", "start", "end")
    .collect()
)

##################################################
# Go through the reads in the bam file
bam = pb.scan_bam(BAM_FILE, tag_fields=["NH"]).select(
    read_id="name",
    source_transcript=pl.col("name").str.split("_").list.get(1).cast(pl.Categorical),
    chrom="chrom",
    start="start",
    end="end",
    cigar="cigar",
    NH="NH",
)

hits = (
    pb.overlap(
        cigar_to_ref_blocks(
            bam,
            other_cols=["read_id", "source_transcript", "NH", "chrom"],
        ),
        exons,
        cols1=["chrom", "block_start", "block_end"],
    )
    .select(
        read_id="read_id_1",
        source_transcript="source_transcript_1",
        multimapper=pl.col("NH_1") > 1,
        mapped_gene_id="gene_id_2",
    )
    .group_by(
        "source_transcript",
        "mapped_gene_id",
        "multimapper",
    )
    .agg(count=pl.col("read_id").n_unique())
)
hits.sink_csv("processed/BEERS2_transcripts_to_STAR_mapped_genes.txt", separator="\t")
