from typing import TYPE_CHECKING

import numpy as np
from sorted_nearest import maketiles, makewindows

from pyranges.names import END_COL, START_COL, TEMP_END_COL, TEMP_START_COL

if TYPE_CHECKING:
    from pyranges import RangeFrame


def _windows(df: "RangeFrame", *, window_size: int, **_) -> "RangeFrame":
    idxs, starts, ends = makewindows(df.index.values, df.Start.values, df.End.values, window_size)

    df = df.reindex(idxs)
    df.loc[:, START_COL] = starts
    df.loc[:, END_COL] = ends

    return df


def _tiles(df: "RangeFrame", tile_size: int, overlap_column: str | None, **_) -> "RangeFrame":
    if overlap_column is not None:
        df = df.copy()
        df.insert(df.shape[1], TEMP_START_COL, df.Start)
        df.insert(df.shape[1], TEMP_END_COL, df.End)

    idxs, starts, ends = maketiles(
        df.index.values,
        df.Start.values,
        df.End.values,
        tile_size,
    )

    df = df.reindex(idxs)
    df.loc[:, START_COL] = starts
    df.loc[:, END_COL] = ends

    if overlap_column is not None:
        overlap = np.minimum(df.End, df[TEMP_END_COL]) - np.maximum(df.Start, df[TEMP_START_COL])
        df.insert(df.shape[1], overlap_column, overlap)
        df = df.drop([TEMP_START_COL, TEMP_END_COL], axis=1)

    return df
