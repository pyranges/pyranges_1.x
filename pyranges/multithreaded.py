from typing import TYPE_CHECKING, Any, Callable, Dict, List, Tuple, Union, Literal

import numpy as np
import pandas as pd
from natsort import natsorted  # type: ignore
from pandas.core.frame import DataFrame

import pyranges as pr
from pyranges.names import (
    CHROM_COL,
    STRAND_COL,
    GENOME_LOC_COLS_WITH_STRAND,
    TEMP_INDEX_COL,
    VALID_STRAND_BEHAVIOR_OPTIONS,
    STRAND_BEHAVIOR_OPPOSITE,
    VALID_STRAND_BEHAVIOR_TYPE,
    STRAND_BEHAVIOR_AUTO,
    STRAND_BEHAVIOR_IGNORE,
)

if TYPE_CHECKING:
    from pyranges.pyranges_main import PyRanges


def merge_dfs(df1: DataFrame, df2: DataFrame) -> DataFrame:
    if not df1.empty and not df2.empty:
        return pd.concat([df1, df2], sort=False).reset_index(drop=True)

    elif df1.empty and df2.empty:
        # can this happen?
        return pd.DataFrame()
    elif df1.empty:
        return df2
    else:
        return df1


def process_results(
    results: List[Any], keys: Union[List[str], List[Tuple[str, str]]]
) -> dict:
    results_dict = {k: r for k, r in zip(keys, results) if r is not None}

    try:
        next(iter(results_dict.values()))
    except StopIteration:  # empty collection
        return results_dict

    # An arbitrary operation might make the keys in the dict and df out of sync.
    # This fixes that by having the PyRanges initializer find the correct keys again..
    try:
        if all(isinstance(v, pd.DataFrame) for v in results_dict.values()):
            df = pd.concat(results_dict.values())
            import pyranges as pr

            _results_dict = pr.PyRanges(df).dfs
        else:
            return results_dict
    except (ValueError, TypeError):
        return results_dict

    to_delete = []
    # to ensure no duplicate indexes and no empty dataframes
    for k in results_dict:
        if results_dict[k] is None or results_dict[k].empty:
            to_delete.append(k)
        else:
            # pandas might make a df that is not always C-contiguous
            # copying fixes this
            # TODO: better to only fix columns that are not C-contiguous?
            results_dict[k] = results_dict[k].copy(deep=True)
            results_dict[k].index = range(len(results_dict[k]))

    for k in to_delete:
        del results_dict[k]

    return _results_dict


def make_sparse(df: DataFrame) -> DataFrame:
    if "Strand" in df:
        cols = "Chromosome Start End Strand".split()
    else:
        cols = "Chromosome Start End".split()

    return df[cols]


def make_binary_sparse(
    kwargs: Dict[str, Any], df: DataFrame, odf: DataFrame
) -> Tuple[DataFrame, DataFrame]:
    sparse = kwargs.get("sparse")

    if not sparse:
        return df, odf

    if sparse.get("self"):
        df = make_sparse(df)

    if sparse.get("other"):
        odf = make_sparse(odf)

    return df, odf


def make_unary_sparse(kwargs: Dict[str, Any], df: DataFrame) -> DataFrame:
    sparse = kwargs.get("sparse", {}).get("self")

    return make_sparse(df) if sparse else df


def _group_keys(
    self: "PyRanges",
    other: "PyRanges",
    strand_behavior: VALID_STRAND_BEHAVIOR_TYPE,
) -> bool:
    include_strand = True
    if strand_behavior == STRAND_BEHAVIOR_AUTO:
        include_strand = self.strand_values_valid and other.strand_values_valid
    elif strand_behavior == STRAND_BEHAVIOR_IGNORE:
        include_strand = False
    return [CHROM_COL, STRAND_COL] if include_strand else [CHROM_COL]


def pyrange_apply(
    function: Callable,
    self: "PyRanges",
    other: "PyRanges",
    strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
    **kwargs,
) -> pd.DataFrame:
    other_strand = {"+": "-", "-": "+"}
    same_strand = {"+": "+", "-": "-"}

    strand_dict = other_strand if strand_behavior == STRAND_BEHAVIOR_OPPOSITE else same_strand

    ensure_strand_options_valid(other, self, strand_behavior)

    results = []

    original_index = self.index.names
    should_reset_index_self = index_contains_genome_loc_cols(self)
    should_reset_index_other = index_contains_genome_loc_cols(other)
    self = self.reset_index() if should_reset_index_self else self
    other = other.reset_index() if should_reset_index_other else other

    grpby_ks = _group_keys(self, other, strand_behavior)

    others = dict(list(other.groupby(grpby_ks, observed=True)))
    empty = pr.empty(columns=other.columns, dtype=other.dtypes)
    for key, gdf in self.groupby(grpby_ks, observed=True):
        if strand_behavior == STRAND_BEHAVIOR_OPPOSITE:
            key = key[0], strand_dict[key[1]]
        ogdf = others.get(key, empty)
        results.append(function(gdf, ogdf, **kwargs))
    result = pd.concat(results)
    return result.set_index(original_index) if should_reset_index_self else result


def ensure_strand_options_valid(other, self, strand_behavior):
    if strand_behavior not in VALID_STRAND_BEHAVIOR_OPTIONS:
        msg = f"{VALID_STRAND_BEHAVIOR_OPTIONS} are the only valid values for strand_behavior."
        raise ValueError(msg)
    if strand_behavior == STRAND_BEHAVIOR_OPPOSITE:
        assert (
            self.strand_values_valid and other.strand_values_valid
        ), "Can only do opposite strand operations when both PyRanges contain valid strand info."


def pyrange_apply_single(
    function: Callable, self: "PyRanges", **kwargs
) -> pd.DataFrame:
    strand = kwargs["strand"]

    if strand and not STRAND_COL in self.columns:
        msg = "Can only do stranded operation when PyRange contains strand info"
        raise ValueError(msg)

    if not self.strand_values_valid:
        keys = [CHROM_COL]
    else:
        keys = [CHROM_COL, STRAND_COL] if strand else [CHROM_COL]
    range_index = np.arange(len(self))
    if isinstance(self.index, pd.RangeIndex):
        self = self.set_index(pd.Series(name=TEMP_INDEX_COL, data=range_index))
        res = (
            self.groupby(keys, as_index=False, observed=True)
            .apply(function, **kwargs)
            .reset_index(drop=True)
        )
        return res
    else:
        if self.index.name is None and self.index.names == [None]:
            original_index = None
        else:
            original_index = (
                self.index.names if self.index.name is None else self.index.names
            )
        self = self.reset_index().set_index(
            pd.Series(name=TEMP_INDEX_COL, data=range_index), append=False
        )
        res = self.groupby(keys, as_index=False, observed=True).apply(
            function, **kwargs
        )
        return (
            res.reset_index(drop=True)
            if original_index is None
            else res.set_index(original_index)
        )


def _lengths(df):
    lengths = df.End - df.Start

    return lengths


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

    assert isinstance(
        slack, (int, dict)
    ), "Extend parameter must be integer or dict, is {}".format(type(slack))

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

    assert (
        df.Start < df.End
    ).all(), "Some intervals are negative or zero length after applying extend!"

    return df


def _extend_grp(df: pd.DataFrame, **kwargs):
    df = df.copy()
    dtype = df.Start.dtype
    slack = kwargs["ext"]
    by = kwargs["group_by"]
    g = df.groupby(by)

    assert isinstance(
        slack, (int, dict)
    ), "Extend parameter must be integer or dict, is {}".format(type(slack))

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

    assert (
        df.Start < df.End
    ).all(), "Some intervals are negative or zero length after applying extend!"

    return df


def index_contains_genome_loc_cols(gr: "PyRanges") -> bool:
    # Reset index if special columns are in the index
    return bool(set(gr.index.names) & set(GENOME_LOC_COLS_WITH_STRAND))
