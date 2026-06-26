import argparse
import polars as pl
import polars_bio as pb
import gzip

parser = argparse.ArgumentParser(
    "strip transcripts from fa if they are not gencode basic and add decoys"
)
parser.add_argument("input", help="input fa.gz")
parser.add_argument("gtf", help="input fa.gz")
parser.add_argument("decoys", help="input fa.gz")
parser.add_argument("output", help="output fa")

args = parser.parse_args()

annot = (
    pb.scan_gtf(args.gtf, attr_fields=["transcript_id", "transcript_version", "tag"])
    .filter(pl.col("type") == "transcript")
    .collect()
)

basic_transcripts = set(
    annot.filter(pl.col("tag").str.contains("basic"))
    .select(
        transcript_id=pl.col("transcript_id")
        + "."
        + pl.col("transcript_version").cast(str)
    )["transcript_id"]
    .unique()
)

pb.sink_fasta(
    pb.scan_fasta(args.input).filter(pl.col("name").is_in(basic_transcripts)),
    args.output,
)

# Now add the decoys on
# A bug in polars_bio prevents us from reading this in with pb.scan_fasta()
# so we just append to the output file
with open(args.output, "at") as out:
    with gzip.open(args.decoys, "rt") as decoys:
        for line in decoys:
            out.write(line)
