from typing import TYPE_CHECKING

import numpy as np

from pyranges.names import (
    END_COL,
    FORWARD_STRAND,
    REVERSE_STRAND,
    START_COL,
    TEMP_END_COL,
    TEMP_INDEX_COL,
    TEMP_MAX_COL,
    TEMP_MIN_COL,
    TEMP_START_COL,
    VALID_GENOMIC_STRAND_TYPE,
)

if TYPE_CHECKING:
    from pyranges import PyRanges


def _subseq(
    scdf: "PyRanges",
    start: int = 0,
    end: int | None = None,
    strand: VALID_GENOMIC_STRAND_TYPE = FORWARD_STRAND,
    **kwargs,
) -> "PyRanges":
    if scdf.empty:
        return None

    scdf = scdf.copy()
    orig_order = scdf.index.copy()
    scdf.insert(0, TEMP_INDEX_COL, orig_order)

    by_argument_given = kwargs.get("by")
    _by = kwargs.get("by", TEMP_INDEX_COL)
    by = [_by] if isinstance(_by, str) else (_by or [])

    # at this point, strand is False if:
    #   1. subsequence was called with strand=False or
    #   2. it was called with strand=None and self is not stranded
    # or it can be '+' or '-' if:
    #   1. it was input as True to subsequence and passed  to pyrange_apply_single as True,
    #      which updates it to '-' or '+' before calling _subseq, or
    #   2. it was called with strand=None and self is stranded
    strand = FORWARD_STRAND if strand != REVERSE_STRAND else strand
    # now, unstranded or strand==None cases are treated like all intervals are on the + strand

    # creating j which holds the boundaries per group
    # j contains one row per group; columns: Start  End (+ by columns); indexed by __i__
    agg_dict = {TEMP_INDEX_COL: "first", START_COL: "min", END_COL: "max"} | {k: "first" for k in by}

    if by_argument_given:
        j = (
            scdf.groupby(by, dropna=False)[[START_COL, END_COL, TEMP_INDEX_COL, *by]]
            .agg(agg_dict)
            .set_index(TEMP_INDEX_COL)
        )
    else:
        j = scdf.groupby(by, dropna=False)[[START_COL, END_COL, TEMP_INDEX_COL]].agg(agg_dict).set_index(TEMP_INDEX_COL)
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
    scdf = scdf.merge(j[[*by, TEMP_MIN_COL, TEMP_MAX_COL]], on=by).set_index(TEMP_INDEX_COL).loc[orig_order]

    # instead of simply using starts and ends as computed above, we're dealing here with potential out of bounds:
    r = scdf[~((scdf[START_COL] >= scdf[TEMP_MAX_COL]) | (scdf[END_COL] <= scdf[TEMP_MIN_COL]))].copy()
    r.loc[:, START_COL] = np.maximum(r.Start, r[TEMP_MIN_COL])
    r.loc[:, END_COL] = np.minimum(r.End, r[TEMP_MAX_COL])

    return r.drop([TEMP_MIN_COL, TEMP_MAX_COL], axis=1)
