import polars as pl

BAM = "/home/thobr/nonuniform_impact/processed/unbiased/samples/S1/bam/Aligned.sortedByCoord.out.bam"

CONSUMES_REF = {
    "M": True,
    "I": False,
    "N": True,
    "D": True,
    "S": False,
    "H": False,
    "P": False,
    "=": True,
    "X": True,
}
INCLUDE_BLOCK = {
    "M": True,
    "I": False,
    "N": False,
    "D": False,  # Do we want to count deletions or not? No, for now but could
    "S": False,
    "H": False,
    "P": False,
    "=": True,
    "X": True,
}


def cigar_to_ref_blocks(
    df, cigar: str = "cigar", start: str = "start", other_cols: list[str] = []
):
    """
    Adds 'block_start' and 'block_end' columns to df and explodes to one row per block in the cigar string

    Computes coordinates of the reference not the query

    EXAMPLE:
        df = pl.DataFrame(
            {"start": [1, 100, 1000], "cigar": ["50M10N50M", "5S10M5S", "100M10N100M"]}
        )
        df.config_meta.set(coordinate_system_zero_based=False)
        cigar_to_ref_blocks(df)
    """
    assert not df.config_meta.get_metadata().get(
        "coordinate_system_zero_based", True
    ), "We only support 1-based coordinates / make sure metadata is set"
    toks = pl.col(cigar).str.extract_all(r"[0-9]+[MINDSHPX=]")
    lens = toks.list.eval(pl.element().str.extract(r"([0-9]+)").cast(pl.Int32))
    ops = toks.list.eval(pl.element().str.extract(r"([MINDSHPX=])"))
    consumes = ops.list.eval(
        pl.element()
        .replace_strict(CONSUMES_REF, return_dtype=pl.Boolean)
        .cast(pl.Int32)
    )
    include = ops.list.eval(
        pl.element().replace_strict(INCLUDE_BLOCK, return_dtype=pl.Boolean)
    )
    local_end = (lens * consumes).list.eval(pl.element().cum_sum())
    local_start = local_end - lens

    return (
        df.with_columns(_ls=local_start, _le=local_end, _inc=include)
        .explode("_ls", "_le", "_inc")
        .filter("_inc")
        .select(
            *other_cols,
            block_start=pl.col(start) + pl.col("_ls"),
            block_end=pl.col(start) + pl.col("_le") - 1,
        )
    )
