import numpy as np
import pandas as pd
from ncls import NCLS  # type: ignore[import]

from pyranges.names import VALID_JOIN_TYPE


def _both_indexes(
    scdf: pd.DataFrame,
    ocdf: pd.DataFrame,
) -> tuple[pd.Index, pd.Index]:
    if ocdf.empty:
        return scdf.index, pd.Index(np.array([], dtype=np.int64))

    starts = scdf.Start.to_numpy()
    ends = scdf.End.to_numpy()
    indexes = scdf.index.to_numpy()

    it = NCLS(ocdf.Start.to_numpy(), ocdf.End.to_numpy(), ocdf.index.to_numpy())

    return it.all_overlaps_both(starts, ends, indexes)


def _both_dfs(
    scdf: pd.DataFrame,
    ocdf: pd.DataFrame,
    join_type: VALID_JOIN_TYPE,
    suffix: str = "_b",
    **_,
) -> pd.DataFrame:
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
        ocdf_missing.index = pd.Index(-np.arange(1, len(ocdf_missing) + 1))
        scdf_missing = scdf.reindex(missing_indices_self)
        scdf_missing.index = pd.Index(-np.arange(len(ocdf_missing) + 1, len(ocdf_missing) + 1 + len(scdf_missing)))
        scdf = pd.concat([scdf_matching, scdf_missing])
        ocdf = pd.concat([ocdf_matcing, ocdf_missing])
    return scdf.merge(
        ocdf,
        left_index=True,
        right_index=True,
        suffixes=("", suffix),
        how="inner" if join_type is None else join_type,
    )
