from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyranges.core.names import (
    CHROM_COL,
    END_COL,
    START_COL,
    STRAND_COL,
)
from pyranges.core.pyranges_helpers import arg_to_list

if TYPE_CHECKING:
    from pyranges import PyRanges

cumsum_start = "_refcumstart"
cumsum_end = "_refcumend"
suffix = "__ref"
start_b = START_COL + suffix
end_b = END_COL + suffix
strand_b = STRAND_COL + suffix


def _map_to_local(gr, ref, idcol, match_by) -> "PyRanges":
    # record ordered columns of gr
    gr_cols = gr.columns.tolist()
    gr_has_strand = gr.has_strand
    ref_has_strand = ref.has_strand
    match_by = arg_to_list(match_by)
    ref_cols_to_keep = [*match_by, idcol]

    if idcol in match_by:
        # dealing with case in which user used idcol as match_by
        fixed_idcol = idcol
        gr_cols = [c for c in gr_cols if c != idcol]
        ref_cols_to_keep.pop(-1)  # remove idcol from match_by
    else:
        fixed_idcol = "__idcol"

    ref = (
        ref.get_with_loc_columns(ref_cols_to_keep)
        .group_cumsum(
            match_by=idcol,
            use_strand="auto",
            cumsum_start_column=cumsum_start,
            cumsum_end_column=cumsum_end,
            sort=False,
        )
        .rename(columns={idcol: fixed_idcol})
    )

    gr = gr.join_ranges(ref, strand_behavior="ignore", match_by=match_by, suffix=suffix)

    # removing gr regions that do not overlap with ref
    gr = gr.combine_interval_columns("intersect", drop_old_columns=False, start2=start_b, end2=end_b)

    columns_to_drop = [cumsum_start, cumsum_end, start_b, end_b, CHROM_COL]

    # dealing with - strand intervals on ref.
    if ref_has_strand:
        # depending on whether gr had strand, STRAND_COL + suffix may be missing in join_ranges output
        strand_col_of_ref = strand_b if gr_has_strand else STRAND_COL

        ref_strand_is_neg = gr[strand_col_of_ref] == "-"

        # computing start and end for neg intervals at the same time, to allow replacement in one go
        gr.loc[ref_strand_is_neg, START_COL], gr.loc[ref_strand_is_neg, END_COL] = (
            (
                gr.loc[ref_strand_is_neg, end_b]
                - gr.loc[ref_strand_is_neg, END_COL]
                + gr.loc[ref_strand_is_neg, cumsum_start]
            ),
            (
                gr.loc[ref_strand_is_neg, end_b]
                - gr.loc[ref_strand_is_neg, START_COL]
                + gr.loc[ref_strand_is_neg, cumsum_start]
            ),
        )

        if gr_has_strand:  # and ref_has_strand: # implicit
            # transform Strand: for rows with Strand == Strand_b, it becomes '+', otherwise '-'
            gr[STRAND_COL] = np.where(gr[STRAND_COL] == gr[strand_b], "+", "-")
            columns_to_drop = [*columns_to_drop, strand_b]
    else:
        ref_strand_is_neg = pd.Series(data=False, index=gr.index, dtype=np.bool_)

    # dealing with + strand intervals in ref
    gr.loc[~ref_strand_is_neg, START_COL] = (
        gr.loc[~ref_strand_is_neg, START_COL]
        - gr.loc[~ref_strand_is_neg, start_b]
        + gr.loc[~ref_strand_is_neg, cumsum_start]
    )
    gr.loc[~ref_strand_is_neg, END_COL] = (
        gr.loc[~ref_strand_is_neg, END_COL]
        - gr.loc[~ref_strand_is_neg, start_b]
        + gr.loc[~ref_strand_is_neg, cumsum_start]
    )

    # drop useless columns
    gr = gr.drop(columns=columns_to_drop).rename(columns={fixed_idcol: CHROM_COL})

    # reordering columns to match original gr
    return gr.reindex(columns=gr_cols)
