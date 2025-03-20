from typing import TYPE_CHECKING

import numpy as np
from ruranges import merge_numpy  # type: ignore[import]

from pyranges.core.names import END_COL, START_COL
from pyranges.core.pyranges_helpers import (
    check_and_return_common_type_2,
    check_min_max_with_slack,
    factorize,
    mypy_ensure_rangeframe,
)

if TYPE_CHECKING:
    from pyranges.range_frame.range_frame import RangeFrame


def _merge(
    df: "RangeFrame",
    by: list[str],
    count_col: str | None = None,
    slack: int | None = None,
) -> "RangeFrame":
    from pyranges.range_frame.range_frame import RangeFrame

    if df.empty:
        return df

    col_order = [col for col in df if col in [*by, START_COL, END_COL]]

    factorized = factorize(df, by)

    starts = df[START_COL].to_numpy()
    ends = df[END_COL].to_numpy()
    _dtype = check_and_return_common_type_2(starts, ends)

    if slack:
        check_min_max_with_slack(starts, ends, slack, _dtype)

    indices, start, end, counts = merge_numpy(
        chrs=factorized,
        starts=df.Start.astype(np.int64).to_numpy(),
        ends=df.End.astype(np.int64).to_numpy(),
        slack=slack,
    )

    by_subset = df[col_order].take(indices)
    by_subset.loc[:, START_COL] = start.astype(_dtype)
    by_subset.loc[:, END_COL] = end.astype(_dtype)

    outpr = RangeFrame(by_subset)

    if count_col:
        outpr.insert(outpr.shape[1], count_col, counts)

    return mypy_ensure_rangeframe(outpr.reset_index(drop=True))
