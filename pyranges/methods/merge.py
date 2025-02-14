from typing import TYPE_CHECKING

from ruranges import merge_numpy  # type: ignore[import]

from pyranges.core.names import END_COL, START_COL
from pyranges.core.pyranges_helpers import factorize, mypy_ensure_rangeframe

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

    indices, start, end, counts = merge_numpy(
        chrs=factorized,
        starts=df.Start.to_numpy(),
        ends=df.End.to_numpy(),
        slack=slack,
    )

    by_subset = df[col_order].take(indices)
    by_subset.loc[:, START_COL] = start
    by_subset.loc[:, END_COL] = end

    outpr = RangeFrame(by_subset)

    if count_col:
        outpr.insert(outpr.shape[1], count_col, counts)

    return mypy_ensure_rangeframe(outpr.reset_index(drop=True))
