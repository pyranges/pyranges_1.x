from typing import TYPE_CHECKING

import pandas as pd
import numpy as np

from pyranges.core.pyranges_helpers import factorize
import ruranges

from pyranges.core.names import (
    END_COL,
    FORWARD_STRAND,
    START_COL,
    STRAND_COL,
)

if TYPE_CHECKING:
    import pyranges as pr



def _spliced_subseq(
    df: "pr.PyRanges",
    *,
    by: list[str],
    force_plus_strand: bool = False,
    start: int = 0,
    end: int | None = None,
) -> pd.DataFrame:
    if df.empty:
        return df

    chrs = factorize(df, by) if by else np.arange(len(df), dtype=np.int64)

    outidx, outstarts, outends = ruranges.spliced_subsequence_numpy(
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
