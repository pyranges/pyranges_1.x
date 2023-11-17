import numpy as np
import pandas as pd
from ncls import NCLS  # type: ignore

from pyranges import empty_df


def _both_indexes(scdf, ocdf, how=False, **kwargs):
    if ocdf.empty:
        return scdf.index, np.array([], dtype=np.int64)

    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    it = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    if not how or how in ["outer", "left", "right"]:
        _self_indexes, _other_indexes = it.all_overlaps_both(starts, ends, indexes)
    elif how == "containment":
        _self_indexes, _other_indexes = it.all_containments_both(starts, ends, indexes)
    elif how == "first":
        _self_indexes, _other_indexes = it.first_overlap_both(starts, ends, indexes)
    elif how == "last":
        _self_indexes, _other_indexes = it.last_overlap_both(starts, ends, indexes)

    return _self_indexes, _other_indexes

def _both_dfs(scdf, ocdf, how=False, **kwargs):
    _self_indexes, _other_indexes = _both_indexes(scdf, ocdf, how, **kwargs)
    if how == "left":
        ocdf = ocdf.reindex(_other_indexes)
        ocdf.index = _self_indexes
    elif how == "right":
        scdf = scdf.reindex(_self_indexes)
        scdf.index = _other_indexes
    elif how == "outer":
        missing_indices_self = scdf.index.difference(_self_indexes)
        missing_indices_other = ocdf.index.difference(_other_indexes)
        scdf_matching = scdf.reindex(_self_indexes)
        scdf_matching.index = _other_indexes
        ocdf_matcing = ocdf.reindex(_other_indexes)
        ocdf_missing = ocdf.reindex(missing_indices_other)
        ocdf_missing.index = -np.arange(1, len(ocdf_missing) + 1)
        scdf_missing = scdf.reindex(missing_indices_self)
        scdf_missing.index = -np.arange(
            len(ocdf_missing) + 1, len(ocdf_missing) + 1 + len(scdf_missing)
        )
        scdf = pd.concat([scdf_matching, scdf_missing])
        ocdf = pd.concat([ocdf_matcing, ocdf_missing])
    return scdf.merge(
        ocdf,
        left_on=scdf.index,
        right_on=ocdf.index,
        suffixes=("", kwargs["suffix"]),
        how="inner" if how is None else how,
    )
