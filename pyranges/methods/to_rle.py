from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import contextlib

    from pyranges import PyRanges

    with contextlib.suppress(ImportError):
        from pyrle import PyRles  # type: ignore[import-not-found]


def _to_rle(
    ranges: "PyRanges",
    value_col: str | None = None,
    *,
    strand: bool = True,
    rpm: bool = False,
    **kwargs,
) -> "PyRles":
    try:
        from pyrle import PyRles, Rle  # type: ignore[import]
        from pyrle.methods import _coverage  # type: ignore[import]
    except ImportError as e:
        msg = "Using the coverage method requires that pyrle is installed."
        raise ImportError(msg) from e

    ranges = ranges.remove_strand() if not strand else ranges
    _kwargs = {
        "strand": strand,
        "value_col": value_col,
        "sparse": {"self": False},
    }  # already sparse
    kwargs.update(_kwargs)

    result = {}
    for k, v in ranges.groupby(ranges.loc_columns):
        value_col = kwargs.get("value_col")
        if value_col:
            values = v[value_col].astype(np.float64).to_numpy(copy=True)
        else:
            values = np.ones(len(v), dtype=np.float64)
        starts = v.Start.to_numpy(copy=True)
        ends = v.End.to_numpy(copy=True)
        positions = np.concatenate([starts, ends])
        weights = np.concatenate([values, -values])
        order = np.argsort(positions, kind="mergesort")
        positions = np.ascontiguousarray(positions[order].astype(np.int64, copy=False))
        weights = np.ascontiguousarray(weights[order])
        runs, outvals = _coverage(positions, weights)
        result[k] = Rle(runs, outvals)

    if rpm:
        multiplier = 1e6 / len(ranges)
        result = {k: v * multiplier for k, v in result.items()}

    return PyRles(result)
