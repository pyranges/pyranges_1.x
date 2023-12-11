from typing import TYPE_CHECKING

import pandas as pd
from ncls import NCLS  # type: ignore[import]

from pyranges.names import TEMP_NUM_COL

if TYPE_CHECKING:
    from pyranges import RangeFrame


def _subtraction(scdf: "RangeFrame", ocdf: "RangeFrame", **_) -> "RangeFrame":
    original_class = scdf.__class__
    if ocdf.empty or scdf.empty:
        return scdf

    o = NCLS(ocdf.Start.to_numpy(), ocdf.End.to_numpy(), ocdf.index.to_numpy())

    idx_self, new_starts, new_ends = o.set_difference_helper(
        scdf.Start.values,
        scdf.End.values,
        scdf.index.values,
        scdf[TEMP_NUM_COL].to_numpy(),
    )

    missing_idx = pd.Index(scdf.index).difference(idx_self)

    idx_to_drop = new_starts != -1

    new_starts = new_starts[idx_to_drop]
    new_ends = new_ends[idx_to_drop]

    idx_self = idx_self[idx_to_drop]
    new_starts = pd.Series(new_starts, index=idx_self)
    new_ends = pd.Series(new_ends, index=idx_self)

    _scdf = scdf.reindex(missing_idx.union(idx_self)).sort_index()
    new_starts = new_starts.sort_index()
    new_ends = new_ends.sort_index()

    if len(idx_self):
        _scdf.loc[_scdf.index.isin(idx_self), "Start"] = new_starts.to_numpy()
        _scdf.loc[_scdf.index.isin(idx_self), "End"] = new_ends.to_numpy()

    return original_class(_scdf)
