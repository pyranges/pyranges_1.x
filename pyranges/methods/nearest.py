from typing import TYPE_CHECKING

import pandas as pd
from ncls import NCLS  # type: ignore[import]
from sorted_nearest import (  # type: ignore[import]
    nearest_next_nonoverlapping,
    nearest_nonoverlapping,
    nearest_previous_nonoverlapping,
)

import pyranges.empty
from pyranges import PyRanges
from pyranges.methods.sort import sort_one_by_one
from pyranges.names import END_COL, START_COL

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray
    from pandas import DataFrame


def _insert_distance(df2: pd.DataFrame, dist: "NDArray | int", suffix: str) -> pd.DataFrame:
    if "Distance" not in df2:
        distance_column_name = "Distance"
    elif "Distance" + suffix not in df2:
        distance_column_name = "Distance" + suffix
    else:
        i = 1
        while "Distance" + str(i) in df2:
            i += 1
        distance_column_name = "Distance" + str(i)

    df2.insert(
        df2.shape[1],
        distance_column_name,
        pd.Series(dist, index=df2.index).fillna(-1).astype(int),
    )

    return df2


def _overlapping_for_nearest(df: pd.DataFrame, df2: pd.DataFrame, suffix: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    nearest_df = pd.DataFrame(columns="Chromosome Start End Strand".split())

    it = NCLS(df2.Start.to_numpy(), df2.End.to_numpy(), df2.index.to_numpy())

    idx_self, idx_other = it.first_overlap_both(
        df[START_COL].to_numpy(),
        df[END_COL].to_numpy(),
        df.index.to_numpy(),
    )
    df2, df22 = df.reindex(idx_self), df2.reindex(idx_other)

    if not df22.empty:
        idxs = df2.index
        original_idx = df.index.copy(deep=True)
        missing_idxs = ~original_idx.isin(idxs)
        missing_overlap = df.index[missing_idxs]

        df_to_find_nearest_in = df.reindex(missing_overlap)

        odf = df2.reindex(df22.index)
        odf.index = idxs
        sdf = df.reindex(idxs)

        nearest_df = sdf.join(odf, rsuffix=suffix)
        nearest_df = _insert_distance(nearest_df, 0, suffix)
    else:
        df_to_find_nearest_in = df

    return nearest_df, df_to_find_nearest_in


def _next_nonoverlapping(
    left_ends: pd.Series,
    right_starts: pd.Series,
    right_indexes: "ArrayLike",
) -> tuple["NDArray", "NDArray"]:
    left_ends = left_ends.sort_values()
    right_starts = right_starts.sort_values()
    r_idx, dist = nearest_next_nonoverlapping(left_ends.to_numpy() - 1, right_starts.to_numpy(), right_indexes)
    r_idx = pd.Series(r_idx, index=left_ends.index).sort_index().to_numpy()
    dist = pd.Series(dist, index=left_ends.index).sort_index().to_numpy()

    return r_idx, dist


def _previous_nonoverlapping(left_starts: pd.Series, right_ends: pd.Series) -> tuple["NDArray", "NDArray"]:
    left_starts = left_starts.sort_values()
    right_ends = right_ends.sort_values()
    r_idx, dist = nearest_previous_nonoverlapping(
        left_starts.to_numpy(),
        right_ends.to_numpy() - 1,
        right_ends.index.to_numpy(),
    )

    r_idx = pd.Series(r_idx, index=left_starts.index).sort_index().to_numpy()
    dist = pd.Series(dist, index=left_starts.index).sort_index().to_numpy()

    return r_idx, dist


def _nearest(df: "DataFrame", df2: "DataFrame", **kwargs) -> pd.DataFrame:
    if df.empty or df2.empty:
        return PyRanges()

    overlap = kwargs["overlap"]
    how = kwargs["how"]
    suffix = kwargs["suffix"]

    if how == "upstream":
        strand = df.Strand.iloc[0]
        how = {"+": "previous", "-": "next"}[strand]
    elif how == "downstream":
        strand = df.Strand.iloc[0]
        how = {"+": "next", "-": "previous"}[strand]

    df2 = df2.reset_index(drop=True)

    if overlap:
        nearest_df, df_to_find_nearest_in = _overlapping_for_nearest(df, df2, suffix)
    else:
        df_to_find_nearest_in = df
        nearest_df = pyranges.empty.empty()

    df = pyranges.empty.empty_df()
    if not df_to_find_nearest_in.empty:
        df_to_find_nearest_in = sort_one_by_one(df_to_find_nearest_in, "Start", "End")
        df2 = sort_one_by_one(df2, "Start", "End")
        df_to_find_nearest_in.index = pd.Index(range(len(df_to_find_nearest_in)))

        if how == "next":
            r_idx, dist = _next_nonoverlapping(df_to_find_nearest_in.End, df2.Start, df2.index.to_numpy())
        elif how == "previous":
            r_idx, dist = _previous_nonoverlapping(df_to_find_nearest_in.Start, df2.End)
        else:
            previous_r_idx, previous_dist = _previous_nonoverlapping(df_to_find_nearest_in.Start, df2.End)

            next_r_idx, next_dist = _next_nonoverlapping(df_to_find_nearest_in.End, df2.Start, df2.index.to_numpy())

            r_idx, dist = nearest_nonoverlapping(previous_r_idx, previous_dist, next_r_idx, next_dist)

        df2 = df2.reindex(r_idx)

        df2.index = df_to_find_nearest_in.index

        df2 = _insert_distance(df2, dist, suffix)

        _r_idx = pd.Series(r_idx, index=df2.index)
        df_to_find_nearest_in = df_to_find_nearest_in.drop(_r_idx[_r_idx == -1].index)

        df = df_to_find_nearest_in.join(df2, rsuffix=suffix)

    if overlap and "df" in locals() and not df.empty and not nearest_df.empty:
        df = pd.concat([nearest_df, df], sort=False)
    elif overlap and not nearest_df.empty:
        df = nearest_df

    return PyRanges(df.drop("Chromosome" + suffix, axis=1))
