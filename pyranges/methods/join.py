import numpy as np
import pandas as pd
from ncls import NCLS  # type: ignore

from pyranges import empty_df
from pyranges.names import VALID_JOIN_TYPE, JOIN_INNER, JOIN_RIGHT


def _both_indexes(
    scdf,
    ocdf,
) -> tuple[np.array, np.array]:
    if ocdf.empty:
        return scdf.index, np.array([], dtype=np.int64)

    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    it = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    return it.all_overlaps_both(starts, ends, indexes)


def _both_dfs(scdf, ocdf, join_type: VALID_JOIN_TYPE, **kwargs):
    _self_indexes, _other_indexes = _both_indexes(scdf, ocdf)
    if join_type == "inner":
        scdf = scdf.reindex(_self_indexes)
        ocdf = ocdf.reindex(_other_indexes)
        ocdf.index = _self_indexes
    elif join_type == "left":
        ocdf = ocdf.reindex(_other_indexes)
        ocdf.index = _self_indexes
    elif join_type == "right":
        scdf = scdf.reindex(_self_indexes)
        scdf.index = _other_indexes
    elif join_type == "outer":
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
    suffixes = ("", kwargs["suffix"]) if kwargs["suffix"] else None
    return scdf.merge(
        ocdf,
        left_index=True,
        right_index=True,
        suffixes=suffixes,
        how="inner" if join_type is None else join_type,
    )
