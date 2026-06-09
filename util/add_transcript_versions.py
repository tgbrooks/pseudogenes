import argparse
import polars as pl
import polars_bio as pb

parser = argparse.ArgumentParser(
    "map transcripts to genes for Salmon (includes transcript id to match index)"
)
parser.add_argument("input", help="gtf file to read (gzipped allowed)")
parser.add_argument("output", help="txt file to write to")

args = parser.parse_args()
gtf = (
    pb.scan_gtf(
        args.input, attr_fields=["gene_id", "transcript_id", "transcript_version"]
    )
    .filter(pl.col("type") == "transcript")
    .collect()
    .select(
        transcript_id=pl.col("transcript_id") + "." + pl.col("transcript_version"),
        gene_id="gene_id",
    )
)

gtf.write_csv(args.output, separator="\t")
