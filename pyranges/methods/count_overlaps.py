from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import ruranges  # type: ignore[import]

from pyranges.core.names import END_COL, START_COL
from pyranges.core.pyranges_helpers import (
    check_min_max_with_slack,
    factorize_binary,
)

if TYPE_CHECKING:
    from pyranges.range_frame.range_frame import RangeFrame


def _count_overlaps(
    df: "RangeFrame",
    df2: "RangeFrame",
    by: list[str],
    slack: int | None = None,
) -> "pd.Series":
    if df.empty:
        return pd.Series()

    f1, f2 = factorize_binary(df, df2, by)

    starts = df[START_COL].to_numpy()
    ends = df[END_COL].to_numpy()

    if slack:
        check_min_max_with_slack(starts, ends, slack, np.int64)

    return ruranges.count_overlaps_numpy(  # type: ignore[attr-defined]
        f1,
        df[START_COL].astype(np.int64).to_numpy(),
        df[END_COL].astype(np.int64).to_numpy(),
        f2,
        df2[START_COL].astype(np.int64).to_numpy(),
        df2[END_COL].astype(np.int64).to_numpy(),
        slack=slack,
    )
