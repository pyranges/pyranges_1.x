import numpy as np
import pandas as pd
from ncls import NCLS  # type: ignore[import]

from pyranges.names import VALID_JOIN_TYPE


def _both_indexes(
    df: pd.DataFrame,
    df2: pd.DataFrame,
) -> tuple[pd.Index, pd.Index]:
    if df2.empty:
        return df.index, pd.Index(np.array([], dtype=np.int64))

    starts = df.Start.to_numpy()
    ends = df.End.to_numpy()
    indexes = df.index.to_numpy()

    it = NCLS(df2.Start.to_numpy(), df2.End.to_numpy(), df2.index.to_numpy())

    return it.all_overlaps_both(starts, ends, indexes)


def _both_dfs(
    df: pd.DataFrame,
    df2: pd.DataFrame,
    join_type: VALID_JOIN_TYPE,
    suffix: str = "_b",
    **_,
) -> pd.DataFrame:
    _self_indexes, _other_indexes = _both_indexes(df, df2)
    if join_type == "inner":
        df = df.reindex(_self_indexes)
        df2 = df2.reindex(_other_indexes)
        df2.index = _self_indexes
    elif join_type == "left":
        df2 = df2.reindex(_other_indexes)
        df2.index = _self_indexes
    elif join_type == "right":
        df = df.reindex(_self_indexes)
        df.index = _other_indexes
    elif join_type == "outer":
        missing_indices_self = df.index.difference(_self_indexes)
        missing_indices_other = df2.index.difference(_other_indexes)
        df_matching = df.reindex(_self_indexes)
        df_matching.index = _other_indexes
        df2_matcing = df2.reindex(_other_indexes)
        df2_missing = df2.reindex(missing_indices_other)
        df2_missing.index = pd.Index(-np.arange(1, len(df2_missing) + 1))
        df_missing = df.reindex(missing_indices_self)
        df_missing.index = pd.Index(-np.arange(len(df2_missing) + 1, len(df2_missing) + 1 + len(df_missing)))
        df = pd.concat([df_matching, df_missing])
        df2 = pd.concat([df2_matcing, df2_missing])
    return df.merge(
        df2,
        left_index=True,
        right_index=True,
        suffixes=("", suffix),
        how="inner" if join_type is None else join_type,
    )
