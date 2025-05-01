from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
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


def map_to_global(  # noqa: C901  (intentionally long)
    local_gr: "PyRanges",
    global_gr: "PyRanges",
    global_on: str,
    *,
    local_on: str = "Chromosome",
    use_strand: VALID_USE_STRAND_TYPE = USE_STRAND_DEFAULT,
) -> "PyRanges":
    """
    Map *local* intervals (e.g. transcript-relative) to absolute genomic
    coordinates, using a Rust kernel for speed.

    ..  (docstring unchanged – abridged here for brevity)  ..
    """

    # ------------------------------------------------------------------
    # 1.  Build transcript-local coordinates for each exon in *global_gr*
    # ------------------------------------------------------------------
    cumsum_start = "_local_start"
    cumsum_end = "_local_end"

    global_cum = global_gr.group_cumsum(
        match_by=global_on,
        use_strand=use_strand,
        cumsum_start_column=cumsum_start,
        cumsum_end_column=cumsum_end,
    ).copy()

    # Preserve true genomic coordinates
    global_cum["__global_start__"] = global_cum[START_COL]
    global_cum["__global_end__"] = global_cum[END_COL]

    # Replace Start/End with transcript-local coord expected by Rust
    global_cum[START_COL] = global_cum[cumsum_start]
    global_cum[END_COL] = global_cum[cumsum_end]

    # ------------------------------------------------------------------
    # 2.  Factorise transcript IDs so exons & queries share *numeric* codes
    # ------------------------------------------------------------------
    local_df = local_gr.copy()
    ex_df = global_cum.copy()

    # one concatenation – same codes array split into two slices
    tx_codes, _ = pd.factorize(
        pd.concat([ex_df[global_on], local_df[local_on]], ignore_index=True),
        sort=False,
    )
    ex_tx_code = tx_codes[: len(ex_df)].astype(np.uint32, copy=False)
    q_tx_code = tx_codes[len(ex_df) :].astype(np.uint32, copy=False)

    # ------------------------------------------------------------------
    # 3.  Factorise genomic *chromosome* names (returned by Rust)
    # ------------------------------------------------------------------
    chr_codes, chr_table = pd.factorize(ex_df[CHROM_COL], sort=False)
    ex_chr_code = chr_codes.astype(np.uint32, copy=False)

    # ------------------------------------------------------------------
    # 4.  Build numpy views for the Rust function
    # ------------------------------------------------------------------
    ex_local_start = ex_df[START_COL].to_numpy(np.int64, copy=False)
    ex_local_end = ex_df[END_COL].to_numpy(np.int64, copy=False)
    ex_genome_start = ex_df["__global_start__"].to_numpy(np.int64, copy=False)
    ex_genome_end = ex_df["__global_end__"].to_numpy(np.int64, copy=False)
    ex_fwd = (ex_df[STRAND_COL] == "+").to_numpy(np.bool_, copy=False)

    q_start = local_df[START_COL].to_numpy(np.int64, copy=False)
    q_end = local_df[END_COL].to_numpy(np.int64, copy=False)
    q_idx = np.arange(len(local_df), dtype=np.uint32)  # 0-based row numbers
    q_fwd = (local_df[STRAND_COL] == "+").to_numpy(np.bool_, copy=False)

    # ------------------------------------------------------------------
    # 5.  Call the compiled Rust kernel
    # ------------------------------------------------------------------
    (
        keep_idx,  # row numbers into *local_df*
        out_chr_code,
        out_start,
        out_end,
        out_strand_bool,
    ) = ruranges.map_to_global_numpy(
        # exon side
        ex_tx_code,
        ex_local_start,
        ex_local_end,
        ex_chr_code,
        ex_genome_start,
        ex_genome_end,
        ex_fwd,
        # query side
        q_tx_code,
        q_start,
        q_end,
        q_fwd,
    )

    # ------------------------------------------------------------------
    # 6.  Nothing mapped?  -> return empty PyRanges with same columns
    # ------------------------------------------------------------------
    if len(keep_idx) == 0:
        return local_gr.__class__(pd.DataFrame(columns=local_df.columns), int64=True)

    # ------------------------------------------------------------------
    # 7.  Re-assemble result by *taking* relevant rows, then patch coords
    # ------------------------------------------------------------------
    mapped = local_gr.take(keep_idx)  # duplicates preserved
    mapped_df = mapped.copy()

    mapped_df[CHROM_COL] = chr_table[out_chr_code]
    mapped_df[START_COL] = out_start
    mapped_df[END_COL] = out_end
    mapped_df[STRAND_COL] = np.where(out_strand_bool, "+", "-")

    return local_gr.__class__(mapped_df)
