from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyranges.core.names import (
    END_COL,
    FORWARD_STRAND,
    START_COL,
    STRAND_COL,
)
from pyranges.core.pyranges_helpers import factorize

if TYPE_CHECKING:
    import pyranges as pr


def _spliced_subseq(
    df: "pr.PyRanges",
    *,
    by: list[str],
    force_plus_strand: bool = False,
    start: int | Sequence[int] | np.ndarray = 0,
    end: int | Sequence[int] | np.ndarray | None = None,
) -> pd.DataFrame:
    import ruranges

    if df.empty:
        return df

    chrs = factorize(df, by) if by else np.arange(len(df), dtype=np.uint32)

    starts = np.ascontiguousarray(df[START_COL].to_numpy())
    ends = np.ascontiguousarray(df[END_COL].to_numpy())
    if STRAND_COL in df and not force_plus_strand:
        strand_flags = np.ascontiguousarray((df[STRAND_COL] == FORWARD_STRAND).to_numpy())
    else:
        strand_flags = np.ones(len(df), dtype=bool)
    if isinstance(start, np.ndarray):
        start = np.ascontiguousarray(start)
    if isinstance(end, np.ndarray):
        end = np.ascontiguousarray(end)

    outidx, outstarts, outends = ruranges.spliced_subsequence(
        groups=chrs,  # type: ignore[arg-type]
        starts=starts,
        ends=ends,
        strand_flags=strand_flags,
        start=start,  # type: ignore[type]
        end=end,  # type: ignore[type]
        force_plus_strand=force_plus_strand,
    )

    outdf = df.take(outidx)  # type: ignore[arg-type]
    outdf.loc[:, START_COL] = outstarts
    outdf.loc[:, END_COL] = outends

    return outdf


# deprecated. Due to unexpected behavior, this is not used anymore (see slice_ranges with count_introns=False)
def _subseq(
    df: "pr.PyRanges",
    *,
    by: list[str],
    start: int = 0,
    end: int | None = None,
    force_plus_strand: bool = False,
) -> pd.DataFrame:
    import ruranges

    if df.empty:
        return df

    chrs = factorize(df, by)

    starts = np.ascontiguousarray(df[START_COL].to_numpy())
    ends = np.ascontiguousarray(df[END_COL].to_numpy())
    if STRAND_COL in df and not force_plus_strand:
        strand_flags = np.ascontiguousarray((df[STRAND_COL] == FORWARD_STRAND).to_numpy())
    else:
        strand_flags = np.ones(len(df), dtype=bool)

    outidx, outstarts, outends = ruranges.subsequence_numpy(  # type: ignore[attr-defined]
        chrs,
        starts,
        ends,
        strand_flags,
        start=start,
        end=end,
        force_plus_strand=force_plus_strand,
    )

    outdf = df.take(outidx)
    outdf.loc[:, START_COL] = outstarts
    outdf.loc[:, END_COL] = outends

    return outdf
