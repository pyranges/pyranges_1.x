import numpy as np
import pandas as pd
from ncls import NCLS  # type: ignore[import]


def _number_overlapping(scdf: pd.DataFrame, ocdf: pd.DataFrame, **kwargs) -> pd.DataFrame:
    keep_nonoverlapping = kwargs.get("keep_nonoverlapping", True)
    column_name = kwargs.get("overlap_col", True)

    if scdf.empty:
        return scdf
    if ocdf.empty:
        if keep_nonoverlapping:
            df = scdf.copy()
            df.insert(df.shape[1], column_name, 0)
            return df
        return ocdf

    oncls = NCLS(ocdf.Start.to_numpy(), ocdf.End.to_numpy(), ocdf.index.to_numpy())

    starts = scdf.Start.to_numpy()
    ends = scdf.End.to_numpy()
    indexes = scdf.index.to_numpy()

    _self_indexes, _other_indexes = oncls.all_overlaps_both(starts, ends, indexes)

    s = pd.Series(_self_indexes)
    counts_per_read = s.value_counts()[s.unique()].reset_index()
    counts_per_read.columns = ["Index", "Count"]

    df = scdf.copy()

    _missing_indexes = np.setdiff1d(scdf.index, _self_indexes)
    missing = pd.DataFrame(data={"Index": _missing_indexes, "Count": 0}, index=_missing_indexes)
    counts_per_read = pd.concat([counts_per_read, missing])

    counts_per_read = counts_per_read.set_index("Index").sort_index().squeeze()

    df.insert(df.shape[1], column_name, counts_per_read)

    if keep_nonoverlapping:
        return df
    return df[df[column_name] != 0]


def _coverage(scdf: pd.DataFrame, ocdf: pd.DataFrame, **kwargs) -> pd.DataFrame:
    fraction_col = kwargs["fraction_col"]

    if scdf.empty:
        return scdf
    if ocdf.empty:
        df = scdf.copy()
        df.insert(df.shape[1], fraction_col, 0.0)
        return df

    oncls = NCLS(ocdf.Start.to_numpy(), ocdf.End.to_numpy(), ocdf.index.to_numpy())

    starts = scdf.Start.to_numpy()
    ends = scdf.End.to_numpy()
    indexes = scdf.index.to_numpy()

    _lengths = oncls.coverage(starts, ends, indexes)
    _lengths = _lengths / (ends - starts)
    _fractions = _lengths
    _fractions = _fractions.astype("float64")
    _fractions = np.nan_to_num(_fractions)

    scdf = scdf.copy()

    scdf.insert(scdf.shape[1], fraction_col, _fractions)

    return scdf
