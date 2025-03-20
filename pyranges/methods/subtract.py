from typing import TYPE_CHECKING

import ruranges

from pyranges.core.names import END_COL, START_COL
from pyranges.core.pyranges_helpers import check_and_return_common_type_4, factorize_binary, mypy_ensure_rangeframe

if TYPE_CHECKING:
    from pyranges.range_frame.range_frame import RangeFrame


def _subtract(
    df: "RangeFrame",
    df2: "RangeFrame",
    *,
    by: list[str],
) -> "RangeFrame":
    f1, f2 = factorize_binary(df, df2, by)

    starts = df[START_COL].to_numpy()
    ends = df[END_COL].to_numpy()
    starts2 = df2[START_COL].to_numpy()
    ends2 = df2[END_COL].to_numpy()
    _dtype = check_and_return_common_type_4(starts, ends, starts2, ends2)

    idx, start, end = ruranges.subtract_numpy(  # type: ignore[attr-defined]
        f1,
        df[START_COL].to_numpy(),
        df[END_COL].to_numpy(),
        f2,
        df2[START_COL].to_numpy(),
        df2[END_COL].to_numpy(),
    )

    output = df.take(idx).copy()
    output[START_COL], output[END_COL] = start.astype(_dtype), end.astype(_dtype)
    return mypy_ensure_rangeframe(output)
