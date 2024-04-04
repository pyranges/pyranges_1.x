from typing import TYPE_CHECKING

import numpy as np
from sorted_nearest import maketiles, makewindows  # type: ignore[import-untyped]

from pyranges.names import BY_ENTRY_IN_KWARGS, END_COL, START_COL, TEMP_END_COL, TEMP_START_COL

if TYPE_CHECKING:
    from pyranges import RangeFrame


def _windows(
    df: "RangeFrame",
    *,
    window_size: int,
    **kwargs,
) -> "RangeFrame":
    strand = kwargs.get(BY_ENTRY_IN_KWARGS, {}).get("Strand")
    # at this point, strand is False/None if 1. function was called with strand=False or
    #                                   2. it was called with strand=None and self is not stranded
    # or it can be '+' or '-' if:
    #  1. it was input as True to function and passed  to pyrange_apply_single as True,
    #     which updates it to '-' or '+' before calling _windows, or
    #  2. it was called with strand=None and self is stranded

    # for neg strand, we use the same makewindows function, so we need to reverse the coordinates
    if strand == "-":
        orig_starts = df.End.to_numpy() * -1
        orig_ends = df.Start.to_numpy() * -1
    else:
        orig_starts = df.Start.to_numpy()
        orig_ends = df.End.to_numpy()

    idxs, starts, ends = makewindows(
        df.index.to_numpy(),
        orig_starts,
        orig_ends,
        window_size,
    )

    if strand == "-":
        starts, ends = (
            ends * -1,
            starts * -1,
        )

    _df = df.reindex(idxs)

    _df.loc[:, START_COL] = starts
    _df.loc[:, END_COL] = ends

    return _df


def _tiles(
    df: "RangeFrame",
    *,
    tile_size: int,
    overlap_column: str | None,
    **_,
) -> "RangeFrame":
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

    _df = df.reindex(idxs)
    _df.loc[:, START_COL] = starts
    _df.loc[:, END_COL] = ends

    if overlap_column is not None:
        overlap = np.minimum(_df.End, _df[TEMP_END_COL]) - np.maximum(_df.Start, _df[TEMP_START_COL])
        _df.insert(_df.shape[1], overlap_column, overlap)
        _df = _df.drop_and_return([TEMP_START_COL, TEMP_END_COL], axis=1)

    return _df
