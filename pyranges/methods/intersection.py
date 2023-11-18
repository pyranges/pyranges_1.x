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


def _intersection(scdf, ocdf, **kwargs):
    how = kwargs["how"]

    if ocdf.empty or scdf.empty:
        return None

    assert how in "containment first last".split() + [False, None]
    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    in_dtype = ocdf.Start.dtype

    oncls = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    if not how or how is None:
        _self_indexes, _other_indexes = oncls.all_overlaps_both(starts, ends, indexes)
    elif how == "containment":
        _self_indexes, _other_indexes = oncls.all_containments_both(
            starts, ends, indexes
        )
    elif how == "first":
        _self_indexes, _other_indexes = oncls.first_overlap_both(starts, ends, indexes)
    elif how == "last":
        _self_indexes, _other_indexes = oncls.last_overlap_both(starts, ends, indexes)

    _self_indexes = _self_indexes
    _other_indexes = _other_indexes

    scdf, ocdf = scdf.reindex(_self_indexes), ocdf.reindex(_other_indexes)

    new_starts = pd.Series(
        np.where(scdf.Start.values > ocdf.Start.values, scdf.Start, ocdf.Start),
        index=scdf.index,
        dtype=in_dtype,
    )

    new_ends = pd.Series(
        np.where(scdf.End.values < ocdf.End.values, scdf.End, ocdf.End),
        index=scdf.index,
        dtype=in_dtype,
    )

    pd.options.mode.chained_assignment = None  # default='warn'
    scdf.loc[:, "Start"] = new_starts
    scdf.loc[:, "End"] = new_ends
    pd.options.mode.chained_assignment = "warn"

    if not scdf.empty:
        return scdf
    else:
        return None


def _overlap_indices(
    scdf, ocdf, how: Literal["first", "containment", "all"] = "first", **_
) -> "pd.Series[np.int64]":
    if scdf.empty or ocdf.empty:
        return None

    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    it = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    if how == OVERLAP_ALL:
        _indexes = it.all_overlaps_self(starts, ends, indexes)
    elif how == OVERLAP_CONTAINMENT:
        _indexes, _ = it.all_containments_both(starts, ends, indexes)
    elif how == OVERLAP_FIRST:
        _indexes = it.has_overlaps(starts, ends, indexes)
    else:
        msg = f"{VALID_OVERLAP_OPTIONS} are the only valid values for how."
        raise ValueError(msg)
    return _indexes


def _overlap(
    scdf: "pr.RangeFrame",
    ocdf: "pr.RangeFrame",
    how: Literal["first", "containment", "all"] = "first",
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
        return _result.drop("__ix__")
    return result


def _count_overlaps(scdf, ocdf, **kwargs) -> "pd.Series[int]":
    # FIXME: Modifies in-place.
    kwargs["return_indexes"] = True
    idx = _overlap(scdf, ocdf, **kwargs)

    sx = pd.DataFrame(np.zeros(len(scdf), dtype=np.int32), index=scdf.index)

    return pd.Series(idx, dtype=np.int32).value_counts(sort=False)
