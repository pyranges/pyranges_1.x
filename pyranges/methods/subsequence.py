from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyranges.core.names import (
    BY_ENTRY_IN_KWARGS,
    CHROM_COL,
    END_COL,
    FORWARD_STRAND,
    REVERSE_STRAND,
    START_COL,
    STRAND_COL,
    TEMP_END_COL,
    TEMP_INDEX_COL,
    TEMP_MAX_COL,
    TEMP_MIN_COL,
    TEMP_START_COL,
)

if TYPE_CHECKING:
    from pyranges import PyRanges


def _subseq(
    df: "PyRanges",
    *,
    by: list[str],
    start: int = 0,
    end: int | None = None,
    **kwargs,
) -> pd.DataFrame:
    if df.empty:
        return df

    df = df.copy()
    orig_order = df.index.copy()
    df.insert(0, TEMP_INDEX_COL, orig_order)

    by = [col for col in by if col not in [CHROM_COL, STRAND_COL]]
    by_argument_given = bool(by)
    by = by or [TEMP_INDEX_COL]

    strand = kwargs.get(BY_ENTRY_IN_KWARGS, {}).get("Strand")

    # at this point, strand is False if 1. spliced_subsequence was called with use_strand=False or
    #                                   2. it was called with use_strand='auto' and not self.valid_strand
    # or it can be '+' or '-' if:
    #  1. it was input as True to spliced_subsequence and passed  to pyrange_apply_single as True,
    #     which updates it to '-' or '+' before calling _spliced_subseq, or
    #  2. it was called with strand=None and self is stranded

    strand = strand if strand else FORWARD_STRAND

    # creating j which holds the boundaries per group
    # j contains one row per group; columns: Start  End (+ by columns); indexed by __i__
    agg_dict = {TEMP_INDEX_COL: "first", START_COL: "min", END_COL: "max"} | {k: "first" for k in by}

    if by_argument_given:
        j = (
            df.groupby(by, dropna=False)[[START_COL, END_COL, TEMP_INDEX_COL, *by]]
            .agg(agg_dict)
            .set_index(TEMP_INDEX_COL)
        )
    else:
        j = df.groupby(by, dropna=False)[[START_COL, END_COL, TEMP_INDEX_COL]].agg(agg_dict).set_index(TEMP_INDEX_COL)
        j.insert(0, TEMP_INDEX_COL, j.index)
        j.index.name = None

    _end = end if end is not None else (j.End - j.Start).max()

    # below we add columns starts, ends to j
    # start and ends define the desired start and end, with one row per group (transcript).
    # they may be out of bounds of the interval, though (this is dealt with later)

    if (strand == FORWARD_STRAND and start >= 0) or (strand == REVERSE_STRAND and start < 0):
        j.loc[:, TEMP_START_COL] = j.Start + abs(start)
    else:
        j.loc[:, TEMP_START_COL] = j.End - abs(start)

    if (strand == FORWARD_STRAND and _end >= 0) or (strand == REVERSE_STRAND and _end < 0):
        j.loc[:, TEMP_END_COL] = j.Start + abs(_end)
    else:
        j.loc[:, TEMP_END_COL] = j.End - abs(_end)

    if strand == REVERSE_STRAND:
        j = j.rename(columns={TEMP_END_COL: TEMP_MIN_COL, TEMP_START_COL: TEMP_MAX_COL})
    else:
        j = j.rename(columns={TEMP_START_COL: TEMP_MIN_COL, TEMP_END_COL: TEMP_MAX_COL})

    # I'm maintaing the original row order
    _df = df.merge(j[[*by, TEMP_MIN_COL, TEMP_MAX_COL]], on=by).set_index(TEMP_INDEX_COL).loc[orig_order]

    # instead of simply using starts and ends as computed above, we're dealing here with potential out of bounds:
    r = _df[~((_df[START_COL] >= _df[TEMP_MAX_COL]) | (_df[END_COL] <= _df[TEMP_MIN_COL]))].copy()
    r.loc[:, START_COL] = np.maximum(r.Start, r[TEMP_MIN_COL])
    r.loc[:, END_COL] = np.minimum(r.End, r[TEMP_MAX_COL])

    return r.drop([TEMP_MIN_COL, TEMP_MAX_COL], axis=1)
