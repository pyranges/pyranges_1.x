from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from ncls import NCLS  # type: ignore[import]

from pyranges.names import (
    OVERLAP_ALL,
    OVERLAP_CONTAINMENT,
    OVERLAP_FIRST,
    VALID_OVERLAP_OPTIONS,
)

if TYPE_CHECKING:
    import pyranges as pr
    from pyranges import RangeFrame


def _overlap_indices(
    scdf: "RangeFrame",
    ocdf: "RangeFrame",
    how: Literal["first", "containment", "all"] = "first",
    **_,
) -> "pd.Series[np.int64]":
    if scdf.empty or ocdf.empty:
        return np.array([], dtype=np.int64)

    starts = scdf.Start.to_numpy()
    ends = scdf.End.to_numpy()
    indexes = scdf.index.to_numpy()

    it = NCLS(ocdf.Start.to_numpy(), ocdf.End.to_numpy(), ocdf.index.to_numpy())

    if how == OVERLAP_ALL:
        _indices = it.all_overlaps_self(starts, ends, indexes)
    elif how == OVERLAP_CONTAINMENT:
        _indices, _ = it.all_containments_both(starts, ends, indexes)
    elif how == OVERLAP_FIRST:
        _indices, _2 = it.first_overlap_both(starts, ends, indexes)
    else:
        msg = f"{VALID_OVERLAP_OPTIONS} are the only valid to_numpy() for how."
        raise ValueError(msg)
    return _indices


def _overlap(
    scdf: "pr.RangeFrame",
    ocdf: "pr.RangeFrame",
    how: VALID_OVERLAP_OPTIONS = "all",
    *,
    invert: bool = False,
    return_indexes: bool = False,
) -> "pr.RangeFrame":
    if invert:
        scdf = scdf.copy()
        scdf.insert(scdf.shape[1], "__ix__", np.arange(len(scdf)))

    indexes = _overlap_indices(scdf, ocdf, how)
    if return_indexes:
        return indexes
    result = scdf.reindex(indexes)

    if invert:
        found_idxs = getattr(result, "__ix__", [])
        _result = scdf[~pd.Series(scdf.__ix__).isin(found_idxs)]
        return _result.drop("__ix__", axis=1)
    return result


def _count_overlaps(scdf: "RangeFrame", ocdf: "RangeFrame", name: str, return_indexes: bool, **_) -> "pd.Series[int]":
    # FIXME: Modifies in-place.
    idx = _overlap(scdf, ocdf, return_indexes=return_indexes)
    vc = pd.Series(idx, dtype=np.int32).value_counts(sort=False)
    sx = pd.DataFrame(np.zeros(len(scdf), dtype=np.int64), index=scdf.index)

    sx.loc[vc.index, 0] = vc.to_numpy()

    scdf.insert(scdf.shape[1], name, sx)
    return scdf
