from typing import TYPE_CHECKING

import pandas as pd

from pyranges.names import CHROM_COL, END_COL, START_COL, STRAND_COL

if TYPE_CHECKING:
    import pyranges as pr


def _bounds[T: ("pr.PyRanges", "pd.DataFrame")](df: T, **kwargs) -> pd.DataFrame:
    if df.empty:
        return df

    col_order = list(df.columns)

    by = kwargs.get("by")
    by = [by] if isinstance(by, str) else (by or [])

    agg_dict = agg if (agg := kwargs.get("agg")) else {}
    agg_dict.update({START_COL: "min", END_COL: "max", CHROM_COL: "first"})
    if STRAND_COL in df.columns:
        agg_dict[STRAND_COL] = "first"

    res = df.groupby(by).agg(agg_dict).reset_index()
    return res.reindex(columns=[c for c in col_order if c in res.columns])
