from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import ruranges

from pyranges.core.names import (
    END_COL,
    FORWARD_STRAND,
    START_COL,
    STRAND_COL,
)
from pyranges.core.pyranges_helpers import check_and_return_common_type_2, factorize

if TYPE_CHECKING:
    from pyranges import PyRanges


def _subseq(
    df: "PyRanges",
    *,
    by: list[str],
    start: int = 0,
    end: int | None = None,
    force_plus_strand: bool = False,
) -> pd.DataFrame:
    if df.empty:
        return df

    starts = df[START_COL].to_numpy()
    ends = df[END_COL].to_numpy()
    _dtype = check_and_return_common_type_2(starts, ends)

    factorized = factorize(df, by)

    outidx, outstarts, outends = ruranges.subsequence_numpy(  # type: ignore[attr-defined]
        factorized,
        starts.astype(np.int64),
        ends.astype(np.int64),
        (df[STRAND_COL] == FORWARD_STRAND).to_numpy()
        if (STRAND_COL in df and not force_plus_strand)
        else np.ones(len(df), dtype=bool),
        start=start,
        end=end,
        force_plus_strand=force_plus_strand,
    )

    outdf = df.take(outidx)
    outdf.loc[:, START_COL] = outstarts.astype(_dtype)
    outdf.loc[:, END_COL] = outends.astype(_dtype)

    return outdf
