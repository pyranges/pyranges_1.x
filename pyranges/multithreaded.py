from collections.abc import Callable
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame

from pyranges.names import (
    END_COL,
    START_COL,
    STRAND_COL,
    VALID_BY_TYPES,
    VALID_STRAND_TYPE,
)
from pyranges.pyranges_helpers import group_keys_single

if TYPE_CHECKING:
    from pyranges.pyranges_main import PyRanges


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
    **kwargs,
) -> DataFrame:
    df = df.copy()
    dtype = df.Start.dtype

    if ext is not None:
        df.loc[:, "Start"] = df.Start - ext
        df.loc[df.Start < 0, "Start"] = 0
        df.End = df.End + ext
    else:
        if len(strands := df.Strand.drop_duplicates()) > 1:
            msg = f"Cannot extend intervals with different strands: {strands}"
            raise ValueError(msg)
        strand = strands.iloc[0]

        if ext_5 and strand == "+":
            df.loc[:, "Start"] -= ext_5
        elif ext_5 and strand == "-":
            df.loc[:, "End"] += ext_5

        if ext_3 and strand == "-":
            df.loc[:, "Start"] -= ext_3
        elif ext_3 and strand == "+":
            df.loc[:, "End"] += ext_3

    df = df.astype({"Start": dtype, "End": dtype})

    if not (df.Start < df.End).all():
        msg = "Some intervals are negative or zero length after applying extend!"
        raise ValueError(msg)

    return df


def _extend_grp(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    df = df.copy()
    dtype = df.Start.dtype
    slack = kwargs["ext"]
    by = kwargs["group_by"]
    g = df.groupby(by)

    if not isinstance(slack, int | dict):
        msg = f"Extend parameter must be integer or dict, is {type(slack)}"
        raise TypeError(msg)

    minstarts_pos = g.Start.idxmin()
    maxends_pos = g.End.idxmax()

    if isinstance(slack, int):
        df.loc[minstarts_pos, START_COL] = df.Start - slack
        df.loc[df.Start < 0, START_COL] = 0
        df.loc[maxends_pos, END_COL] = df.End + slack

    else:
        strand = df.Strand.iloc[0]
        five_end_slack = slack.get("5")
        three_end_slack = slack.get("3")

        if five_end_slack and strand == "+":
            df.loc[minstarts_pos, START_COL] -= five_end_slack
        elif five_end_slack and strand == "-":
            df.loc[maxends_pos, END_COL] += five_end_slack

        if three_end_slack and strand == "-":
            df.loc[minstarts_pos, START_COL] -= three_end_slack
        elif three_end_slack and strand == "+":
            df.loc[maxends_pos, END_COL] += three_end_slack

    df = df.astype({START_COL: dtype, END_COL: dtype})

    if not (df.Start < df.End).all():
        msg = "Some intervals are negative or zero length after applying extend!"
        raise ValueError(msg)

    return df
