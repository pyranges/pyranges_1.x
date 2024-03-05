from typing import TYPE_CHECKING

import pandas as pd

from pyranges.names import (
    BY_ENTRY_IN_KWARGS,
    END_COL,
    REVERSE_STRAND,
    START_COL,
    TEMP_CUMSUM_COL,
    TEMP_INDEX_COL,
    TEMP_LENGTH_COL,
)

if TYPE_CHECKING:
    import pyranges as pr


TOTAL_EXON_LENGTH_COL = "__temp_total_exon_len__"


def _spliced_subseq(
    df: "pr.PyRanges",
    *,
    start: int = 0,
    end: int | None = None,
    **kwargs,
) -> pd.DataFrame:
    if df.empty:
        return df

    df = df.copy()

    by_argument_given = kwargs.get("by")
    _by = kwargs.get("by", TEMP_INDEX_COL)
    by = [_by] if isinstance(_by, str) else (_by or df.index)

    strand = kwargs.get(BY_ENTRY_IN_KWARGS, {}).get("Strand")

    # at this point, strand is False if 1. spliced_subsequence was called with strand=False or
    #                                   2. it was called with strand=None and self is not stranded
    # or it can be '+' or '-' if:
    #  1. it was input as True to spliced_subsequence and passed  to pyrange_apply_single as True,
    #     which updates it to '-' or '+' before calling _spliced_subseq, or
    #  2. it was called with strand=None and self is stranded

    if strand and not df.strand_values_valid:
        msg = "Cannot have strand=True on unstranded pyranges!"
        raise AssertionError(msg)

    df.insert(df.shape[1], TEMP_LENGTH_COL, df.End - df.Start)
    df.insert(df.shape[1], TEMP_INDEX_COL, df.index)

    g = df.groupby(by, dropna=False)
    df.insert(df.shape[1], TEMP_CUMSUM_COL, g[TEMP_LENGTH_COL].cumsum())

    end = df[TEMP_CUMSUM_COL].max() if end is None else end

    minstart_idx = g[TEMP_INDEX_COL].first()

    if start < 0 or (end is not None and end < 0):
        # len_per_transc is total sum of exon length per transcript
        len_per_transc = df.loc[g[TEMP_INDEX_COL].last(), [*by, TEMP_CUMSUM_COL]].rename(
            columns={TEMP_CUMSUM_COL: TOTAL_EXON_LENGTH_COL},
        )

        # exp_len_per_transc has same rows of df with total sum of exon length
        # had to add bits to keep the order of rows right, or merge would destroy it
        if by_argument_given:
            exp_len_per_transc = (
                df.loc[:, [*by, TEMP_INDEX_COL]].merge(len_per_transc, on=by).set_index(TEMP_INDEX_COL).loc[df.index]
            )
        else:
            exp_len_per_transc = df.loc[:, by].merge(len_per_transc, on=by).set_index(TEMP_INDEX_COL).loc[df.index]  # type: ignore[reportArgumentType, reportCallIssue]

        starts = exp_len_per_transc[TOTAL_EXON_LENGTH_COL] + start if start < 0 else start

        ends = exp_len_per_transc[TOTAL_EXON_LENGTH_COL] + end if end < 0 else end
    else:
        starts = start
        ends = end

    cs_start = g[TEMP_CUMSUM_COL].shift(1, fill_value=0)
    cs_start.loc[minstart_idx] = 0

    cs_end = df[TEMP_CUMSUM_COL]

    # NOTE
    # here below, start is a scalar if originally provided > 0, or a Series if < 0
    #             end is a scalar if originally None or provided >0, or a Series if provided < 0
    if strand == REVERSE_STRAND:  # and use_strand:
        start_adjustments = starts - cs_start
        adjust_start = start_adjustments > 0
        df.loc[adjust_start, END_COL] -= start_adjustments[adjust_start].astype(df[END_COL].dtype)

        end_adjustments = cs_end - ends
        adjust_end = end_adjustments > 0
        df.loc[adjust_end, START_COL] += end_adjustments[adjust_end].astype(df[START_COL].dtype)
    else:
        start_adjustments = starts - cs_start
        adjust_start = start_adjustments > 0
        df.loc[adjust_start, START_COL] += start_adjustments[adjust_start].astype(df[START_COL].dtype)

        end_adjustments = cs_end - ends
        adjust_end = end_adjustments > 0
        df.loc[adjust_end, END_COL] -= end_adjustments[adjust_end].astype(df[END_COL].dtype)

    return df[(df[START_COL] < df[END_COL])].drop([TEMP_INDEX_COL, TEMP_LENGTH_COL, TEMP_CUMSUM_COL], axis=1)
