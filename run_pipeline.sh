#!/usr/bin/env sh
set -e
module load /project/itmatlab/sharedmodules/STAR-v2.7.10b
module load ncbi-blast/2.14.0
module load samtools/1.20

bsub -e logs/snakemake.err \
    -o logs/snakemake.out \
    uv run snakemake --executor lsf \
    --default-resources mem_mb=4000 lsf_queue=normal \
    -j 100 -c 100 \
    --use-singularity --singularity-args "-B /project/itmatlab/" \
    "$@"
