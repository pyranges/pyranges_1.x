from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from ruranges import complement_numpy  # type: ignore[import]

from pyranges.core.names import CHROM_COL, END_COL, START_COL
from pyranges.core.pyranges_helpers import (
    check_and_return_common_type_2,
    check_min_max_with_slack,
    factorize,
    mypy_ensure_rangeframe,
)

if TYPE_CHECKING:
    from pyranges.range_frame.range_frame import RangeFrame


def _complement(
    df: "RangeFrame",
    *,
    by: list[str],
    slack: int = 0,
    chromsizes_col: str | None = None,
    chromsizes: "dict[str | int, int] | None" = None,
    include_first_interval: bool = False,
) -> "RangeFrame":
    from pyranges.range_frame.range_frame import RangeFrame

    if df.empty:
        return df

    starts = df[START_COL].to_numpy()
    ends = df[END_COL].to_numpy()
    _dtype = check_and_return_common_type_2(starts, ends)

    col_order = [col for col in df if col in [*by, START_COL, END_COL]]

    factorized = factorize(df, by)

    if slack:
        check_min_max_with_slack(starts, ends, slack, _dtype)

    if chromsizes and chromsizes_col:
        chrom_lens = df[chromsizes_col].replace(chromsizes)
        chrom_lens = pd.concat([pd.Series(factorized), chrom_lens], axis=1).drop_duplicates()
        chrom_len_ids = chrom_lens[0].to_numpy()
        chrom_lens = chrom_lens[chromsizes_col].to_numpy()
    else:
        chrom_len_ids = np.array([], dtype=np.uint32)
        chrom_lens = np.array([], dtype=np.int64)

    chrs, start, end, idxs = complement_numpy(
        chrs=factorized,
        starts=df.Start.astype(np.int64).to_numpy(),
        ends=df.End.astype(np.int64).to_numpy(),
        slack=slack,
        chrom_len_ids=chrom_len_ids,
        chrom_lens=chrom_lens,
        include_first_interval=include_first_interval,
    )

    ids = df.take(idxs)

    result = RangeFrame(
        {CHROM_COL: chrs, START_COL: start.astype(_dtype), END_COL: end.astype(_dtype)} | {_by: ids[_by] for _by in by}
    )[col_order]

    return mypy_ensure_rangeframe(result.reset_index(drop=True))
