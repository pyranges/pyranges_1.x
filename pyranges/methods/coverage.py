import numpy as np
import pandas as pd
from ncls import NCLS  # type: ignore[import]


def _number_overlapping(df: pd.DataFrame, df2: pd.DataFrame, **kwargs) -> pd.DataFrame:
    keep_nonoverlapping = kwargs.get("keep_nonoverlapping", True)
    column_name = kwargs.get("overlap_col", True)

    if df.empty:
        return df
    if df2.empty:
        if keep_nonoverlapping:
            df = df.copy()
            df.insert(df.shape[1], column_name, 0)
            return df
        return df2

    oncls = NCLS(df2.Start.to_numpy(), df2.End.to_numpy(), df2.index.to_numpy())

    starts = df.Start.to_numpy()
    ends = df.End.to_numpy()
    indexes = df.index.to_numpy()

    _self_indexes, _other_indexes = oncls.all_overlaps_both(starts, ends, indexes)

    s = pd.Series(_self_indexes)
    counts_per_read = s.value_counts()[s.unique()].reset_index()
    counts_per_read.columns = ["Index", "Count"]

    df = df.copy()

    _missing_indexes = np.setdiff1d(df.index, _self_indexes)
    missing = pd.DataFrame(data={"Index": _missing_indexes, "Count": 0}, index=_missing_indexes)
    counts_per_read = pd.concat([counts_per_read, missing])

    counts_per_read = counts_per_read.set_index("Index").sort_index().squeeze()

    df.insert(df.shape[1], column_name, counts_per_read)

    if keep_nonoverlapping:
        return df
    return df[df[column_name] != 0]


def _coverage(df: pd.DataFrame, df2: pd.DataFrame, **kwargs) -> pd.DataFrame:
    fraction_col = kwargs["fraction_col"]

    if df.empty:
        return df
    if df2.empty:
        df = df.copy()
        df.insert(df.shape[1], fraction_col, 0.0)
        return df

    oncls = NCLS(df2.Start.to_numpy(), df2.End.to_numpy(), df2.index.to_numpy())

    starts = df.Start.to_numpy()
    ends = df.End.to_numpy()
    indexes = df.index.to_numpy()

    _lengths = oncls.coverage(starts, ends, indexes)
    _lengths = _lengths / (ends - starts)
    _fractions = _lengths
    _fractions = _fractions.astype("float64")
    _fractions = np.nan_to_num(_fractions)

    df = df.copy()

    df.insert(df.shape[1], fraction_col, _fractions)

    return df
