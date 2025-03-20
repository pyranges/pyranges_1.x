from typing import TYPE_CHECKING

import numpy as np
import ruranges

from pyranges.core.names import END_COL, START_COL
from pyranges.core.pyranges_helpers import check_min_max_with_slack, factorize, mypy_ensure_rangeframe

if TYPE_CHECKING:
    from pyranges import RangeFrame


def _cluster(
    df: "RangeFrame",
    by: list[str],
    cluster_column: str | None = None,
    slack: int | None = None,
) -> "RangeFrame":
    starts = df[START_COL].to_numpy()
    ends = df[END_COL].to_numpy()

    factorized = factorize(df, by)

    if slack:
        check_min_max_with_slack(starts, ends, slack, np.int64)

    cluster, idx = ruranges.cluster_numpy(  # type: ignore[attr-defined]
        factorized,
        starts.astype(np.int64),
        ends.astype(np.int64),
        slack,
    )

    res = df.take(idx).copy()
    res.insert(res.shape[1], cluster_column, cluster)
    return mypy_ensure_rangeframe(res)
