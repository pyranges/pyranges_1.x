from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from pyranges.core.pyranges_helpers import factorize_binary
import ruranges  # Rust extension

from pyranges.core.names import (
    CHROM_COL,
    END_COL,
    START_COL,
    STRAND_COL,
    USE_STRAND_DEFAULT,
    VALID_USE_STRAND_TYPE,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyranges import PyRanges


def map_to_global(
    local_gr: "PyRanges",
    global_gr: "PyRanges",
    global_on: str,
    *,
    local_on: str = "Chromosome",
    use_strand: VALID_USE_STRAND_TYPE = USE_STRAND_DEFAULT,
) -> "PyRanges":
    """
    Lift intervals in *local_gr* (coordinates relative to a transcript)
    onto absolute genomic coordinates using the exon annotation in
    *global_gr*.

    ...
    """

    # 1) Build transcript-local exon coordinates
    cumsum_start = "_local_start"
    cumsum_end   = "_local_end"

    global_cum = global_gr.group_cumsum(
        match_by=global_on,
        use_strand=use_strand,
        cumsum_start_column=cumsum_start,
        cumsum_end_column=cumsum_end,
    )
    ex_df = global_cum.copy()  # DataFrame-like

    ex_df["__global_start__"] = ex_df[START_COL]
    ex_df["__global_end__"]   = ex_df[END_COL]
    ex_df[START_COL] = ex_df[cumsum_start]
    ex_df[END_COL]   = ex_df[cumsum_end]

    # 2) Factorize the transcript IDs in both dataframes
    TMP = "__txid__"
    local_df = local_gr.copy()

    ex_df[TMP]    = ex_df[global_on]
    local_df[TMP] = local_df[local_on]

    q_tx_code, ex_tx_code = factorize_binary(local_df, ex_df, by=TMP)

    # 3) Factorize the chromosome strings in `ex_df`
    ex_chr_code, chr_table = pd.factorize(ex_df[CHROM_COL], sort=False)
    ex_chr_code = ex_chr_code.astype(ex_tx_code.dtype, copy=False)

    # 4) Extract raw NumPy arrays
    ex_local_start  = ex_df[START_COL].to_numpy(copy=False)
    ex_local_end    = ex_df[END_COL].to_numpy(copy=False)
    ex_genome_start = ex_df["__global_start__"].to_numpy(copy=False)
    ex_genome_end   = ex_df["__global_end__"].to_numpy(copy=False)
    ex_fwd          = (ex_df[STRAND_COL] == "+").to_numpy(np.bool_, copy=False)

    q_start = local_df[START_COL].to_numpy(copy=False)
    q_end   = local_df[END_COL].to_numpy(copy=False)
    q_fwd   = (local_df[STRAND_COL] == "+").to_numpy(np.bool_, copy=False)

    # 5) Call the generic Rust kernel
    import ruranges

    # *KEY CHANGE* here:  we put the “exons” as the LEFT triple,
    # and the “queries” as the RIGHT triple to match the Rust function.
    keep_idx, out_start, out_end, out_strand_bool = ruranges.map_to_global(
        # exons => left
        groups=ex_tx_code,
        starts=ex_local_start,
        ends=ex_local_end,
        strand=ex_fwd,
        # queries => right
        groups2=q_tx_code,
        starts2=q_start,
        ends2=q_end,
        chr_code2=ex_chr_code,
        genome_start2=ex_genome_start,
        genome_end2=ex_genome_end,
        strand2=q_fwd,
    )
    # (Alternatively, if you prefer: left=exons, right=queries is consistent
    #  with typical "_dispatch_binary" usage.)

    if len(keep_idx) == 0:
        return pr.PyRanges(local_df.iloc[:0].copy())

    # 6) Convert code → chromosome
    max_code = int(ex_tx_code.max()) + 1
    tx2chr_code = np.empty(max_code, dtype=ex_chr_code.dtype)
    tx2chr_code[ex_tx_code] = ex_chr_code
    out_chr_code = tx2chr_code[q_tx_code[keep_idx]]
    out_chr_str  = chr_table[out_chr_code]

    # 7) Assemble final PyRanges
    mapped = local_gr.take(keep_idx)
    mapped_df = mapped.copy()

    mapped_df[CHROM_COL] = out_chr_str
    mapped_df[START_COL] = out_start
    mapped_df[END_COL]   = out_end
    mapped_df[STRAND_COL] = np.where(out_strand_bool, "+", "-")

    if TMP in mapped_df.columns:
        mapped_df = mapped_df.drop(columns=[TMP])

    return mapped_df