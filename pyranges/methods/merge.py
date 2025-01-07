from pathlib import Path
from typing import TYPE_CHECKING
import pandas as pd
from pyranges.core.pyranges_helpers import mypy_ensure_pyranges
from ruranges import merge_numpy  # type: ignore[import]

from pyranges.core.names import CHROM_COL, END_COL, START_COL

if TYPE_CHECKING:
    import pyranges as pr


def _merge(
    df: "pr.PyRanges",
    by: list[str],
    count_col: str | None = None,
    slack: int | None = None,
) -> "pr.PyRanges":
    from pyranges.core.pyranges_main import PyRanges

    if df.empty:
        return df

    col_order = [col for col in df if col in by + [START_COL, END_COL]]

    # _slack = slack or 0
    chrs = pd.DataFrame(df).groupby(by).ngroup()

    indices, start, end, counts = merge_numpy(
        chrs=chrs.values,
        starts=df.Start.to_numpy(),
        ends=df.End.to_numpy(),
        idxs=df.index.to_numpy(),
        slack=slack,
    )

    by_subset = df[col_order].loc[indices]
    by_subset.loc[:, START_COL] = start
    by_subset.loc[:, END_COL] = end

    outpr = PyRanges(by_subset)

    if count_col:
        outpr.insert(outpr.shape[1], count_col, counts)

    return mypy_ensure_pyranges(outpr.reset_index(drop=True))
