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
    suffix: str,
    **_,
) -> pd.DataFrame:
    _self_indexes, _other_indexes = _both_indexes(df, df2)
    expected_columns = [*df.head(0).join(df2.head(0), how="inner", rsuffix=suffix).columns]
    df2 = pd.DataFrame(df2)
    df2.columns = expected_columns[df.shape[1] :]
    _df = df.reindex(_self_indexes)
    _df.index = pd.Index(np.arange(len(_df)))
    _df2 = pd.DataFrame(df2).reindex(_other_indexes)
    _df2.index = pd.Index(np.arange(len(_df2)))
    j = _df.join(_df2, how="inner")
    if join_type == "inner":
        j.index = _self_indexes
        return j

    if join_type == "left":
        missing_rows = _missing_rows_left(_self_indexes, df, j)
        return pd.concat([j, missing_rows])

    if join_type == "right":
        missing_rows = _missing_rows_right(_other_indexes, df2, j)
        return pd.concat([j, missing_rows])

    if join_type == "outer":
        missing_rows_l = _missing_rows_left(_self_indexes, df, j)
        missing_rows_r = _missing_rows_right(_other_indexes, df2, j)
        return pd.concat([j, missing_rows_l, missing_rows_r])

    msg = f"Invalid join type: {join_type}"
    raise ValueError(msg)


def _missing_rows_left(_self_indexes, df, j) -> pd.DataFrame:
    missing_indices_left = df.index.difference(_self_indexes)
    j.index = _self_indexes
    return df.reindex(missing_indices_left)


def _missing_rows_right(_other_indexes, df2, j) -> pd.DataFrame:
    missing_indices_right = df2.index.difference(_other_indexes)
    j.index = _other_indexes
    return df2.reindex(missing_indices_right)
