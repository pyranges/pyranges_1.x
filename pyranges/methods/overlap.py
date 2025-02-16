from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import ruranges
from numpy.typing import NDArray


from pyranges.core.names import (
    END_COL,
    OVERLAP_ALL,
    OVERLAP_CONTAINED,
    OVERLAP_FIRST,
    OVERLAP_LAST,
    RANGE_COLS,
    START_COL,
    VALID_BY_TYPES,
    VALID_OVERLAP_OPTIONS,
    VALID_OVERLAP_TYPE,
)
from pyranges.core.pyranges_helpers import factorize_binary

if TYPE_CHECKING:
    import pyranges as pr
    from pyranges import RangeFrame


def _overlap_indices(
    df: "RangeFrame",
    df2: "RangeFrame",
    how: VALID_OVERLAP_TYPE = "first",
    **_,
) -> "pd.Series[int]":
    if df.empty or df2.empty:
        return pd.Series([], dtype=np.int64)

    starts = df.Start.to_numpy()
    ends = df.End.to_numpy()
    indexes = df.index.to_numpy()

    it = NCLS(df2.Start.to_numpy(), df2.End.to_numpy(), df2.index.to_numpy())

    if how == OVERLAP_ALL:
        _indices = it.all_overlaps_self(starts, ends, indexes)
    elif how == OVERLAP_CONTAINMENT:
        _indices, _ = it.all_containments_both(starts, ends, indexes)
    elif how == OVERLAP_FIRST:
        _indices, _ = it.first_overlap_both(starts, ends, indexes)
    elif how == OVERLAP_LAST:
        _indices, _ = it.last_overlap_both(starts, ends, indexes)
    else:
        msg = f"{VALID_OVERLAP_OPTIONS} are the only valid values for how."
        raise ValueError(msg)
    return pd.Series(_indices)


def _overlapping_indices(
    df: "pr.RangeFrame",
    df2: "pr.RangeFrame",
    how: VALID_OVERLAP_TYPE = "all",
) -> "pd.Series[int]":
    return _overlap_indices(df, df2, how)


def _both_idxs(
    df: "RangeFrame",
    df2: "RangeFrame",
    *,
    by: VALID_BY_TYPES,
    multiple: VALID_OVERLAP_TYPE = "all",
    contained: bool = False,
    slack: int = 0,
) -> tuple[NDArray[np.int_], NDArray[np.int_]]:
    f1, f2 = factorize_binary(df, df2, by)

    idx1, idx2 = ruranges.chromsweep_numpy(
        f1,
        df.Start.values,
        df.End.values,
        f2,
        df2.Start.values,
        df2.End.values,
        slack,
        overlap_type=multiple,
        contained=contained,
    )
    return idx1, idx2


def _overlap(
    df: "RangeFrame",
    df2: "RangeFrame",
    *,
    by: VALID_BY_TYPES,
    multiple: VALID_OVERLAP_TYPE = "all",
    contained: bool = False,
    slack: int = 0,
) -> pd.DataFrame:
    idx1, _ = _both_idxs(
        df=df,
        df2=df2,
        by=by,
        multiple=multiple,
        contained=contained,
        slack=slack,
    )
    result = df.take(idx1)

    return result


def _count_overlaps(
    df: "RangeFrame",
    df2: "RangeFrame",
    name: str,
    **_,
) -> pd.DataFrame:
    idx = _overlapping_indices(df, df2)
    vc: "pd.Series[int]" = idx.value_counts(sort=False)
    sx = pd.DataFrame(np.zeros(len(df), dtype=np.int64), index=df.index)

    sx.loc[vc.index, 0] = vc.to_numpy()

    df.insert(df.shape[1], name, sx.squeeze())
    return df


def _intersect(
    df: "RangeFrame",
    df2: "RangeFrame",
    *,
    by: list[str],
    multiple: VALID_OVERLAP_TYPE = "all",
    slack: int = 0,
    contained: bool = False,
) -> pd.DataFrame:
    import ruranges

    f1, f2 = factorize_binary(df, df2, by)

    idx1, idx2 = ruranges.chromsweep_numpy(
        f1,
        df.Start.values,
        df.End.values,
        f2,
        df2.Start.values,
        df2.End.values,
        slack,
        overlap_type=multiple,
        contained=contained,
    )

    rf, rf2 = df.take(idx1), df2.take(idx2).loc[:, RANGE_COLS]

    new_starts = np.where(rf.Start.to_numpy() > rf2.Start.to_numpy(), rf.Start, rf2.Start)

    new_ends = np.where(rf.End.to_numpy() < rf2.End.to_numpy(), rf.End, rf2.End)

    rf.loc[:, START_COL] = new_starts
    rf.loc[:, END_COL] = new_ends

    return rf
