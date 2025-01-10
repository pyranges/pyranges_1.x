from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyranges.core.pyranges_helpers import factorize
import ruranges

from pyranges.core.names import (
    END_COL,
    FORWARD_STRAND,
    START_COL,
    STRAND_COL,
)

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

    chrs = factorize(df, by)

    outidx, outstarts, outends = ruranges.subsequence_numpy(
        chrs,
        df[START_COL].to_numpy(),
        df[END_COL].to_numpy(),
        df.index.to_numpy(),
        (df[STRAND_COL] == FORWARD_STRAND).to_numpy() if (STRAND_COL in df and not force_plus_strand) else np.ones(len(df), dtype=bool),
        start=start,
        end=end,
        force_plus_strand=force_plus_strand,
    )

    outdf = df.loc[outidx]
    outdf.loc[:, START_COL] = outstarts
    outdf.loc[:, END_COL] = outends

    return outdf
