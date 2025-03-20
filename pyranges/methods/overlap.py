from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import ruranges
from numpy.typing import NDArray

from pyranges.core.names import (
    END_COL,
    RANGE_COLS,
    START_COL,
    VALID_BY_TYPES,
    VALID_OVERLAP_TYPE,
)
from pyranges.core.pyranges_helpers import check_and_return_common_type_4, check_min_max_with_slack, factorize_binary

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
    f1, f2 = factorize_binary(df, df2, by)
    starts = df[START_COL].to_numpy()
    ends = df[END_COL].to_numpy()
    starts2 = df2[START_COL].to_numpy()
    ends2 = df2[END_COL].to_numpy()
    _dtype = check_and_return_common_type_4(starts, ends, starts2, ends2)

    if slack:
        check_min_max_with_slack(starts, ends, slack, _dtype)

    idx1, idx2 = ruranges.chromsweep_numpy(  # type: ignore[attr-defined]
        f1.astype(np.uint32),
        starts.astype(np.int64),
        ends.astype(np.int64),
        f2.astype(np.uint32),
        starts2.astype(np.int64),
        ends2.astype(np.int64),
        slack,
        overlap_type=multiple,
        contained=contained,
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
        df,
        df2,
        by=by,
        multiple=multiple,
        contained=contained,
        slack=slack,
    )

    rf, rf2 = df.take(idx1), df2.take(idx2).loc[:, RANGE_COLS]

    new_starts = np.where(rf.Start.to_numpy() > rf2.Start.to_numpy(), rf.Start, rf2.Start)

    new_ends = np.where(rf.End.to_numpy() < rf2.End.to_numpy(), rf.End, rf2.End)

    rf.loc[:, START_COL] = new_starts
    rf.loc[:, END_COL] = new_ends

    return rf
