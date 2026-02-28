from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import contextlib

    from pyranges1 import PyRanges

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
        from pyrle import PyRles  # type: ignore[import]
        from pyrle.methods import Rle, _coverage  # type: ignore[import]
    except ImportError as e:
        msg = "Using the coverage method requires that pyrle is installed."
        raise ImportError(msg) from e

    def _coverage_with_writable_arrays(df: pd.DataFrame, *, value_col: str | None = None) -> "Rle":
        values = df[value_col].astype(np.float64).to_numpy(copy=True) if value_col else np.ones(len(df))

        starts_df = pd.DataFrame({"Position": df.Start, "Value": values})[["Position", "Value"]]
        ends_df = pd.DataFrame({"Position": df.End, "Value": -1 * values})[["Position", "Value"]]
        coverage_df = pd.concat([starts_df, ends_df], ignore_index=True)
        coverage_df = coverage_df.sort_values("Position", kind="mergesort")
        coverage_df.loc[:, "Position"] = coverage_df.Position.astype(np.int64)

        # pandas 3 may expose read-only arrays due copy-on-write.
        runs, values = _coverage(
            coverage_df.Position.to_numpy(copy=True),
            coverage_df.Value.to_numpy(copy=True),
        )
        return Rle(runs, values)

    ranges = ranges.remove_strand() if not strand else ranges
    _kwargs = {
        "strand": strand,
        "value_col": value_col,
        "sparse": {"self": False},
    }  # already sparse
    kwargs.update(_kwargs)

    result = {
        k: _coverage_with_writable_arrays(v, value_col=kwargs.get("value_col"))
        for k, v in ranges.groupby(ranges.loc_columns)
    }

    if rpm:
        multiplier = 1e6 / len(ranges)
        result = {k: v * multiplier for k, v in result.items()}

    return PyRles(result)
