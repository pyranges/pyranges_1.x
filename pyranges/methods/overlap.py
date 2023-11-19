from typing import Literal

import numpy as np
import pandas as pd
from ncls import NCLS  # type: ignore

from pyranges.names import (
    VALID_OVERLAP_OPTIONS,
    OVERLAP_FIRST,
    OVERLAP_ALL,
    OVERLAP_CONTAINMENT,
)


def _overlap_indices(
    scdf,
    ocdf,
    how: Literal["first", "containment", "all"] = "first",
    **_,
) -> "pd.Series[np.int64]":
    if scdf.empty or ocdf.empty:
        return np.array([], dtype=np.int64)

    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    it = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    if how == OVERLAP_ALL:
        _indices = it.all_overlaps_self(starts, ends, indexes)
    elif how == OVERLAP_CONTAINMENT:
        _indices, _ = it.all_containments_both(starts, ends, indexes)
    elif how == OVERLAP_FIRST:
        _indices, _2 = it.first_overlap_both(starts, ends, indexes)
    else:
        msg = f"{VALID_OVERLAP_OPTIONS} are the only valid values for how."
        raise ValueError(msg)
    return _indices


def _overlap(
    scdf: "pr.RangeFrame",
    ocdf: "pr.RangeFrame",
    how: VALID_OVERLAP_OPTIONS = "first",
    *,
    invert: bool = False,
    **_,
) -> "pr.RangeFrame":
    if invert:
        scdf = scdf.copy()
        scdf.insert(scdf.shape[1], "__ix__", np.arange(len(scdf)))

    result = scdf.reindex(_overlap_indices(scdf, ocdf, how))

    if invert:
        found_idxs = getattr(result, "__ix__", [])
        _result = scdf[~pd.Series(scdf.__ix__).isin(found_idxs)]
        return _result.drop("__ix__", axis=1)
    return result


def _count_overlaps(scdf, ocdf, **kwargs) -> "pd.Series[int]":
    # FIXME: Modifies in-place.
    kwargs["return_indexes"] = True
    idx = _overlap(scdf, ocdf, **kwargs)

    return pd.Series(idx, dtype=np.int32).value_counts(sort=False)
