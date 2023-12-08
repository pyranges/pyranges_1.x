from typing import TYPE_CHECKING

from pyranges.names import CHROM_COL, END_COL, START_COL, STRAND_COL

if TYPE_CHECKING:
    import pyranges as pr


def _bounds(scdf: "pr.PyRanges", **kwargs) -> "pr.PyRanges":
    if scdf.empty:
        return None

    col_order = [c for c in scdf.columns]

    by = kwargs.get("by")
    by = [by] if isinstance(by, str) else (by or [])

    agg_dict = agg if (agg := kwargs.get("agg")) else {}
    agg_dict.update({START_COL: "min", END_COL: "max", CHROM_COL: "first"})
    if STRAND_COL in scdf.columns:
        agg_dict[STRAND_COL] = "first"

    res = scdf.groupby(by).agg(agg_dict).reset_index()
    return res.reindex(columns=[c for c in col_order if c in res.columns])
