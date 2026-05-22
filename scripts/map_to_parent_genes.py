import json
import gzip
import subprocess

import polars as pl
import simple_gtf

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


biotype_summary = (
    annot.filter(feature="gene")["gene_biotype"]
    .value_counts()
    .select(
        "gene_biotype",
        gene_count="count",
    )
    .with_columns(
        pseudogene=pl.col("gene_biotype").cast(str).str.contains("pseudogene")
        # this is an old pseudogene classification listed by ENSEMBL but doesn't appear to be used
        | (pl.col("gene_biotype") == "retrotransposed")
    )
    .sort("gene_count")
)

biotype_summary.write_csv("data/biotype_summary.txt", separator="\t")

## Read in the transcriptome
transcriptome_sequences = {}
for fasta_file in config["HUMAN_TRANSCRIPTOME_FA"]:
    with gzip.open(fasta_file, "rt") as fasta:
        working_id = None
        working_seq = []
        for line in fasta:
            if line.startswith(">"):
                if working_id:
                    transcriptome_sequences[working_id] = "".join(working_seq)
                working_id = line.removeprefix(">").split()[0]
                working_seq = []
            else:
                working_seq.append(line.strip())
        if working_id:
            transcriptome_sequences[working_id] = "".join(working_seq)
            working_id = None

#### Go through the pseudogenes and attempt to align them to determine which gene they come from
# Generate a fastq file of segments from exons of psedogenes
pseudogene_biotypes = sorted(
    biotype_summary.filter("pseudogene")["gene_biotype"].unique()
)
pseudogenes = (
    annot.filter(feature="gene")
    .select("gene_id", "gene_biotype")
    .filter(pl.col("gene_biotype").is_in(pseudogene_biotypes))
)
pseudogene_transcripts = annot.filter(
    pl.col("feature") == "transcript",
    pl.col("gene_id").is_in(list(pseudogenes["gene_id"])),
).with_columns(
    full_tx_id=pl.col("transcript_id") + "." + pl.col("transcript_version").cast(str)
)
pseudogene_exons = annot.filter(
    pl.col("feature") == "exon",
    pl.col("gene_id").is_in(list(pseudogenes["gene_id"])),
)

# Check overlap of annotation and sequences to verify we have what we want
assert pseudogene_transcripts.select(
    have_all_sequences=pl.col("full_tx_id").is_in(transcriptome_sequences.keys()).all()
)["have_all_sequences"][0]

PSEUDOGENE_SEQ_FILE = "processed/pseudogene_sequence.fa.gz"
with gzip.open(PSEUDOGENE_SEQ_FILE, "wt") as pseudogene_seq_file:
    for tx_id in pseudogene_transcripts["full_tx_id"]:
        seq = transcriptome_sequences[tx_id]
        pseudogene_seq_file.write(f">{tx_id}\n")
        pseudogene_seq_file.write(seq)
        pseudogene_seq_file.write("\n")

######################################################################
# Align by blastn to the full transcriptome to get similarity to genes
BLASTN_RESULTS_FILE = "processed/pseudogene_blastn.csv"
blastn_results = []
for i, tx_file in enumerate(config["HUMAN_TRANSCRIPTOME_FA"]):
    res = subprocess.run(
        f'''bash -c "blastn -subject <(zcat {tx_file}) -query <(zcat {PSEUDOGENE_SEQ_FILE}) -outfmt '10 qseqid sseqid pident nident gaps evalue qseq sseq sstrand qstart qend sstart send length'"''',
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
    )
    blastn_results.append(res.stdout.decode())
with open(BLASTN_RESULTS_FILE, "wt") as out:
    out.write(
        "qseqid,sseqid,pident,nident,gaps,evalue,qseq,sseq,sstrand,qstart,qend,sstart,send,length\n"
    )
    for res in blastn_results:
        out.write(res)

#############
# Parse blast results
blastn_results = (
    pl.read_csv(BLASTN_RESULTS_FILE)
    .with_columns(
        pseudogene_id=pl.col("qseqid").str.split(".").list.get(0),
        transcript_id=pl.col("sseqid").str.split(".").list.get(0),
    )
    .join(
        annot.filter(feature="transcript").select(
            transcript_id="transcript_id",
            gene_id="gene_id",
            biotype="transcript_biotype",
        ),
        on="transcript_id",
        how="left",
    )
    .with_columns(
        biotype_is_pseudogene=pl.col("biotype").cast(str).str.contains("pseudogene")
    )
)

##############################
# Select the best gene matches
# Prioritize protein coding, then lncRNA and takes the best evalue
best_matches = (
    blastn_results.filter(~pl.col("biotype_is_pseudogene"))
    .sort(
        "pseudogene_id",
        pl.col("biotype") == "protein_coding",
        pl.col("biotype") == "lncRNA",
        -pl.col("evalue"),  # want the lowest evalue
        pl.col("pident"),
        descending=True,
    )
    .group_by("pseudogene_id")
    .agg(
        pl.col("gene_id").first(),
        pl.col("biotype").first(),
        pl.col("evalue").first(),
        pl.col("pident").first(),
        pl.col("nident").first(),
        pl.col("gaps").first(),
        n_gene_matches=pl.col("gene_id").n_unique(),
        all_gene_matches=pl.col("gene_id").unique(),
    )
)

best_matches.with_columns(
    all_gene_matches=pl.col("all_gene_matches").list.join(";"),
).write_csv("results/parent_gene_mapping.txt", separator="\t")


# Alternative approach we're no longer using with STAR
# # Extract out kmers to align
# KMER_FILE = "processed/pseudogene_kmer_file.fa.gz"
# KMER_SIZE = 31
# KMER_SPACE = 5
# with gzip.open(KMER_FILE, "wt") as kmer_file:
#     for tx_id in pseudogene_transcripts["full_tx_id"]:
#         seq = transcriptome_sequences[tx_id]
#         for kmer_start in range(0, len(seq), KMER_SPACE):
#             kmer_stop = kmer_start + KMER_SIZE
#             if kmer_stop > len(seq) + KMER_SPACE:
#                 continue
#             kmer_seq = seq[kmer_start:kmer_stop]
#             kmer_id = f"{tx_id}-{kmer_start}-{kmer_stop}"
#             kmer_file.write(f">{kmer_id}\n")
#             kmer_file.write(kmer_seq)
#             kmer_file.write("\n")
#
# Align with STAR:
# INDEX_DIR = config["HUMAN_STAR_INDEX"]
# TEMPORARILY USE THE OLD INDEX WHILE THE NEW ONE IS BUILDING:
# INDEX_DIR = "/project/itmatlab/index/STAR-2.7.10b_indexes/GRCh38.ensemblv114.151bp/"
# subprocess.run(
#    # NOTE: we up the number of multimapping sites allowed just in case
#    f"STAR --runThreadN 6 --genomeDir {INDEX_DIR} --readFilesIn {KMER_FILE} --readFilesCommand zcat --outFileNamePrefix processed/pseudogene_kmer_alignments/ --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --outSAMunmapped Within --outFilterMultimapNmax 25 --winAnchorMultimapNmax 75",
#    check=True,
#    shell=True,
# )
