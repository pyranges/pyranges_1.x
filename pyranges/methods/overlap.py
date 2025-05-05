from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from pyranges.core.names import (
    END_COL,
    RANGE_COLS,
    START_COL,
    VALID_BY_TYPES,
    VALID_OVERLAP_TYPE,
)
from pyranges.core.pyranges_helpers import factorize_binary

if TYPE_CHECKING:
    from pyranges import RangeFrame


def _both_idxs(
    df: "RangeFrame",
    df2: "RangeFrame",
    *,
    by: VALID_BY_TYPES,
    multiple: VALID_OVERLAP_TYPE = "all",
    contained: bool = False,
    slack: int = 0,
) -> tuple[NDArray[np.int_], NDArray[np.int_]]:
    import ruranges

    f1, f2 = factorize_binary(df, df2, by)

    idx1, idx2 = ruranges.overlaps(
        groups=f1,
        starts=df.Start.to_numpy(),
        ends=df.End.to_numpy(),
        groups2=f2,
        starts2=df2.Start.to_numpy(),
        ends2=df2.End.to_numpy(),
        multiple=multiple,
        contained=contained,
        slack=slack,
    )
    return idx1, idx2


def _overlap(
    df: "RangeFrame",
    df2: "RangeFrame",
    *,
    by: VALID_BY_TYPES,
    multiple: VALID_OVERLAP_TYPE = "all",
    contained: bool = False,
    slack: int = 0,
) -> pd.DataFrame:
    idx1, _ = _both_idxs(
        df=df,
        df2=df2,
        by=by,
        multiple=multiple,
        contained=contained,
        slack=slack,
    )
    return df.take(idx1)  # type: ignore[]


def _intersect(
    df: "RangeFrame",
    df2: "RangeFrame",
    *,
    by: list[str],
    multiple: VALID_OVERLAP_TYPE = "all",
    slack: int = 0,
    contained: bool = False,
) -> pd.DataFrame:
    idx1, idx2 = _both_idxs(
        df=df,
        df2=df2,
        by=by,
        multiple=multiple,
        contained=contained,
        slack=slack,
    )

    rf, rf2 = df.take(idx1), df2.take(idx2).loc[:, RANGE_COLS]  # type: ignore[arg-type]

    new_starts = np.where(rf.Start.to_numpy() > rf2.Start.to_numpy(), rf.Start, rf2.Start)

    new_ends = np.where(rf.End.to_numpy() < rf2.End.to_numpy(), rf.End, rf2.End)

    rf.loc[:, START_COL] = new_starts
    rf.loc[:, END_COL] = new_ends

    return rf
