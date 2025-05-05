from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyranges.core.names import VALID_BY_TYPES, VALID_JOIN_TYPE, VALID_OVERLAP_TYPE
from pyranges.methods.overlap import _both_idxs

if TYPE_CHECKING:
    from pyranges.range_frame.range_frame import RangeFrame


def _both_dfs(
    df: "RangeFrame",
    df2: "RangeFrame",
    *,
    by: VALID_BY_TYPES = None,
    multiple: VALID_OVERLAP_TYPE = "all",
    contained: bool = False,
    slack: int = 0,
    join_type: VALID_JOIN_TYPE,
    suffix: str,
) -> pd.DataFrame:
    _self_indexes, _other_indexes = _both_idxs(
        df=df,
        df2=df2,
        by=by,
        multiple=multiple,
        contained=contained,
        slack=slack,
    )
    expected_columns = [*df.head(0).join(df2.head(0), how="inner", rsuffix=suffix).columns]
    df2.columns = expected_columns[df.shape[1] :]
    _df = df.take(_self_indexes)  # type: ignore[arg-type]
    _df.index = pd.Index(np.arange(len(_df)))
    _df2 = pd.DataFrame(df2).take(_other_indexes)  # type: ignore[arg-type]
    _df2.index = pd.Index(np.arange(len(_df2)))
    j = _df.join(_df2, how="inner")
    if join_type == "inner":
        j.index = pd.Index(_self_indexes)
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
