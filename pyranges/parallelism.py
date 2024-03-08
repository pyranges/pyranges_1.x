import math
from collections.abc import Generator

import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame

from pyranges.names import BY_ENTRY_IN_KWARGS, END_COL, FORWARD_STRAND, REVERSE_STRAND, START_COL


def run_in_parallel(function, dfs: list[DataFrame], nb_cpu: int, *args, **kwargs) -> Generator:
    """Run a function in parallel on a list of DataFrames."""
    # TODO(endbak): Things to consider:  # noqa: FIX002
    #   https://github.com/endrebak/pyranges1_alpha/issues/21
    #   The index will be different here than in the original DataFrame if the function resets the index
    #   Perhaps we should not concat so the function can return anything?
    from joblib import Parallel, delayed

    return Parallel(n_jobs=nb_cpu)(delayed(function)(df, *args, **kwargs) for df in dfs)  # type: ignore[return-type]


def split_df_into_chunks_without_splitting_groups(
    df: DataFrame,
    *,
    by: list[str],
    nb_splits: int,
) -> list[DataFrame]:
    """Split a DataFrame into chunks without splitting groups.

    Args:
    ----
        df: DataFrame to split.
        by: Column(s) to group by.
        nb_splits: Number of splits.

    Returns:
    -------
        List of DataFrames.

    """
    # Calculate the target size for each chunk
    target_chunk_size = math.ceil(df.shape[0] / nb_splits)

    # Assign a group ID and calculate cumsum of rows per group
    df = df.sort_values(by)
    groupby = df.groupby(by)

    chunk_ids = groupby.size().cumsum().floordiv(target_chunk_size)
    group_sizes = groupby.size()
    df = df.assign(__chunk_id__=np.repeat(a=chunk_ids.to_numpy(), repeats=group_sizes.to_numpy()))

    # Split the DataFrame using groupby
    return [group.drop("__chunk_id__", axis="columns") for _, group in df.groupby("__chunk_id__")]


def _lengths(df: DataFrame) -> pd.Series:
    return df[END_COL] - df[START_COL]


def _tss(df: DataFrame, **kwargs) -> DataFrame:
    df = df.copy(deep=True)
    dtype = df.dtypes["Start"]
    slack = kwargs.get("slack", 0)

    starts = np.where(df.Strand == "+", df.Start, df.End - 1)
    ends = starts + slack + 1
    starts = starts - slack
    starts = np.where(starts < 0, 0, starts)

    df.loc[:, "Start"] = starts.astype(dtype)
    df.loc[:, "End"] = ends.astype(dtype)

    return df


def _tes(df: DataFrame, **kwargs) -> DataFrame:
    df = df.copy(deep=True)
    dtype = df.dtypes["Start"]
    slack = kwargs.get("slack", 0)

    starts = np.where(df.Strand == "+", df.End - 1, df.Start)
    ends = starts + 1 + slack
    starts = starts - slack
    starts = np.where(starts < 0, 0, starts)

    df.loc[:, "Start"] = starts.astype(dtype)
    df.loc[:, "End"] = ends.astype(dtype)

    return df


def _extend(
    df: DataFrame,
    ext: int | None = None,
    ext_3: int | None = None,
    ext_5: int | None = None,
    **_,
) -> DataFrame:
    df = df.copy()
    dtype = df.Start.dtype

    if ext is not None:
        df.loc[:, START_COL] = df.Start - ext
        df.loc[df.Start < 0, START_COL] = 0
        df.End = df.End + ext
    else:
        if len(strands := df.Strand.drop_duplicates()) > 1:
            msg = f"Cannot extend intervals with different strands: {strands}"
            raise ValueError(msg)
        strand = strands.iloc[0]

        if ext_5 and strand == "+":
            df.loc[:, START_COL] -= ext_5
        elif ext_5 and strand == "-":
            df.loc[:, END_COL] += ext_5

        if ext_3 and strand == "-":
            df.loc[:, START_COL] -= ext_3
        elif ext_3 and strand == "+":
            df.loc[:, END_COL] += ext_3

    df = df.astype({START_COL: dtype, END_COL: dtype})

    if not (df.Start < df.End).all():
        msg = "Some intervals are negative or zero length after applying extend!"
        raise ValueError(msg)

    return df


def _extend_grp(
    df: DataFrame,
    group_by: list[str] | str,
    ext: int | None = None,
    ext_3: int | None = None,
    ext_5: int | None = None,
    **kwargs,
) -> DataFrame:
    df = df.copy()
    g = df.groupby(group_by)

    minstarts_pos = g.Start.idxmin()
    maxends_pos = g.End.idxmax()

    if ext is not None:
        df.loc[minstarts_pos, START_COL] = df.Start - ext
        df.loc[df.Start < 0, START_COL] = 0
        df.loc[maxends_pos, END_COL] = df.End + ext

    else:
        strand = kwargs.get(BY_ENTRY_IN_KWARGS, {}).get("Strand")

        if ext_5 and strand == FORWARD_STRAND:
            df.loc[minstarts_pos, START_COL] -= ext_5
        elif ext_5 and strand == REVERSE_STRAND:
            df.loc[maxends_pos, END_COL] += ext_5

        if ext_3 and strand == REVERSE_STRAND:
            df.loc[minstarts_pos, START_COL] -= ext_3
        elif ext_3 and strand == FORWARD_STRAND:
            df.loc[maxends_pos, END_COL] += ext_3

    if not (df.Start < df.End).all():
        msg = "Some intervals are negative or zero length after applying extend!"
        raise ValueError(msg)

    return df
