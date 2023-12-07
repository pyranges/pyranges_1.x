from collections.abc import Callable
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame

from pyranges.names import (
    CHROM_COL,
    END_COL,
    GENOME_LOC_COLS_WITH_STRAND,
    START_COL,
    STRAND_COL,
)

if TYPE_CHECKING:
    from pyranges.pyranges_main import PyRanges


def pyrange_apply_single(
    function: Callable,
    self: "PyRanges",
    **kwargs,
) -> pd.DataFrame:
    temp_index_col = "__index_column_for_apply__"
    strand = kwargs["strand"]

    if strand and STRAND_COL not in self.columns:
        msg = "Can only do stranded operation when PyRange contains strand info"
        raise ValueError(msg)

    if not self.strand_values_valid:
        keys = [CHROM_COL]
    else:
        keys = [CHROM_COL, STRAND_COL] if strand else [CHROM_COL]
    range_index = np.arange(len(self))
    if isinstance(self.index, pd.RangeIndex):
        self = self.set_index(pd.Series(name=temp_index_col, data=range_index))
        res = self.groupby(keys, as_index=False, observed=True).apply(function, **kwargs).reset_index(drop=True)
        return res
    else:
        if self.index.name is None and self.index.names == [None]:
            original_index = None
        else:
            original_index = self.index.names if self.index.name is None else self.index.names
        self = self.reset_index().set_index(pd.Series(name=temp_index_col, data=range_index), append=False)
        res = self.groupby(keys, as_index=False, observed=True).apply(function, **kwargs)
        return res.reset_index(drop=True) if original_index is None else res.set_index(original_index)


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


def _extend(df: DataFrame, **kwargs) -> DataFrame:
    df = df.copy()
    dtype = df.Start.dtype
    slack = kwargs["ext"]

    assert isinstance(slack, int | dict), f"Extend parameter must be integer or dict, is {type(slack)}"

    if isinstance(slack, int):
        df.loc[:, "Start"] = df.Start - slack
        df.loc[df.Start < 0, "Start"] = 0
        df.End = df.End + slack
    else:
        strand = df.Strand.iloc[0]
        five_end_slack = slack.get("5")
        three_end_slack = slack.get("3")

        if five_end_slack and strand == "+":
            df.loc[:, "Start"] -= five_end_slack
        elif five_end_slack and strand == "-":
            df.loc[:, "End"] += five_end_slack

        if three_end_slack and strand == "-":
            df.loc[:, "Start"] -= three_end_slack
        elif three_end_slack and strand == "+":
            df.loc[:, "End"] += three_end_slack

    df = df.astype({"Start": dtype, "End": dtype})

    assert (df.Start < df.End).all(), "Some intervals are negative or zero length after applying extend!"

    return df


def _extend_grp(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    df = df.copy()
    dtype = df.Start.dtype
    slack = kwargs["ext"]
    by = kwargs["group_by"]
    g = df.groupby(by)

    assert isinstance(slack, int | dict), f"Extend parameter must be integer or dict, is {type(slack)}"

    minstarts_pos = g.Start.idxmin()
    maxends_pos = g.End.idxmax()

    if isinstance(slack, int):
        df.loc[minstarts_pos, "Start"] = df.Start - slack
        df.loc[df.Start < 0, "Start"] = 0
        df.loc[maxends_pos, "End"] = df.End + slack

    else:
        strand = df.Strand.iloc[0]
        five_end_slack = slack.get("5")
        three_end_slack = slack.get("3")

        if five_end_slack and strand == "+":
            df.loc[minstarts_pos, "Start"] -= five_end_slack
        elif five_end_slack and strand == "-":
            df.loc[maxends_pos, "End"] += five_end_slack

        if three_end_slack and strand == "-":
            df.loc[minstarts_pos, "Start"] -= three_end_slack
        elif three_end_slack and strand == "+":
            df.loc[maxends_pos, "End"] += three_end_slack

    df = df.astype({"Start": dtype, "End": dtype})

    assert (df.Start < df.End).all(), "Some intervals are negative or zero length after applying extend!"

    return df


def index_contains_genome_loc_cols(gr: "PyRanges") -> bool:
    # Reset index if special columns are in the index
    return bool(set(gr.index.names) & set(GENOME_LOC_COLS_WITH_STRAND))
