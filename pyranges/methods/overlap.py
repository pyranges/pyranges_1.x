from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import ruranges
from ncls import NCLS  # type: ignore[import]

from pyranges.core.names import (
    END_COL,
    OVERLAP_ALL,
    OVERLAP_FIRST,
    OVERLAP_LAST,
    START_COL,
    VALID_OVERLAP_OPTIONS,
    VALID_OVERLAP_TYPE,
)
from pyranges.core.pyranges_helpers import factorize_multiple

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


def _overlap(
    df: "RangeFrame",
    df2: "RangeFrame",
    *,
    by: list[str],
    multiple: VALID_OVERLAP_TYPE = "all",
    contained: bool = False,
    slack: int = 0,
) -> pd.DataFrame:
    f1, f2 = factorize_multiple(df, df2, by)

    idx1, idx2 = ruranges.chromsweep_numpy(
        f1,
        df.Start.values,
        df.End.values,
        df.index.values,
        f2,
        df2.Start.values,
        df2.End.values,
        df2.index.values,
        slack,
    )

    idxs = pd.Series(idx1, index=idx2)
    if multiple in ["first", "last"]:
        idxs = idxs.drop_duplicates(keep=multiple)

    if contained:
        result = df.loc[idxs.to_numpy()]
        df2_pos = df2.loc[idxs.index.to_numpy(), [START_COL, END_COL]]
        result = result[
            (result[START_COL].to_numpy() >= df2_pos[START_COL].to_numpy())
            & (result[END_COL].to_numpy() <= df2_pos[END_COL].to_numpy())
        ]
    else:
        result = df.loc[idxs.to_numpy()]

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
    how: VALID_OVERLAP_TYPE = "all",
    **_,
) -> pd.DataFrame:
    df = df.copy()
    starts = df.Start.to_numpy()
    ends = df.End.to_numpy()
    indexes = df.index.to_numpy()

    in_dtype = df2.Start.dtype

    oncls = NCLS(df2.Start.to_numpy(), df2.End.to_numpy(), df2.index.to_numpy())

    if how == OVERLAP_ALL:
        _self_indexes, _other_indexes = oncls.all_overlaps_both(starts, ends, indexes)
    elif how == OVERLAP_CONTAINMENT:
        _self_indexes, _other_indexes = oncls.all_containments_both(starts, ends, indexes)
    elif how == OVERLAP_FIRST:
        _self_indexes, _other_indexes = oncls.first_overlap_both(starts, ends, indexes)
    elif how == OVERLAP_LAST:
        _self_indexes, _other_indexes = oncls.last_overlap_both(starts, ends, indexes)
    else:
        msg = "Invalid overlap type. Valid types are: first, containment, all, last."
        raise ValueError(msg)

    df, df2 = df.reindex(_self_indexes), df2.reindex(_other_indexes)

    new_starts = pd.Series(
        np.where(df.Start.to_numpy() > df2.Start.to_numpy(), df.Start, df2.Start),
        index=df.index,
        dtype=in_dtype,
    )

    new_ends = pd.Series(
        np.where(df.End.to_numpy() < df2.End.to_numpy(), df.End, df2.End),
        index=df.index,
        dtype=in_dtype,
    )

    df.loc[:, START_COL] = new_starts
    df.loc[:, END_COL] = new_ends

    return df
