import numpy as np
from sorted_nearest import maketiles, makewindows  # type: ignore

from pyranges.names import TEMP_START_COL, TEMP_END_COL, END_COL, START_COL


def _windows(df, *, window_size: int, **_):
    idxs, starts, ends = makewindows(
        df.index.values, df.Start.values, df.End.values, window_size
    )

    df = df.reindex(idxs)
    df.loc[:, START_COL] = starts
    df.loc[:, END_COL] = ends

    return df


def _tiles(df: "pr.RangeFrame", tile_size: int, overlap_column: str | None, **_):

    if overlap_column is not None:
        df = df.copy()
        df.insert(df.shape[1], TEMP_START_COL, df.Start)
        df.insert(df.shape[1], TEMP_END_COL, df.End)

    idxs, starts, ends = maketiles(
        df.index.values, df.Start.values, df.End.values, tile_size,
    )

    df = df.reindex(idxs)
    df.loc[:, START_COL] = starts
    df.loc[:, END_COL] = ends

    if overlap_column is not None:
        overlap = np.minimum(df.End, df[TEMP_END_COL]) - np.maximum(df.Start, df[TEMP_START_COL])
        df.insert(df.shape[1], overlap_column, overlap)
        df = df.drop([TEMP_START_COL, TEMP_END_COL], axis=1)

    return df
