from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyranges.core.names import CHROM_COL, END_COL, START_COL
from pyranges.core.pyranges_helpers import check_and_return_common_type_2, factorize, mypy_ensure_rangeframe
from pyranges.range_frame.range_frame import RangeFrame

if TYPE_CHECKING:
    import pyranges as pr


def _bounds[T: ("pr.PyRanges", "pd.DataFrame")](df: T, by: list[str]) -> pd.DataFrame:
    if df.empty:
        return df

    col_order = [*dict.fromkeys([col for col in df if col in [*by, START_COL, END_COL]]).keys()]

    import ruranges

    group_ids = factorize(df, by=by)

    starts = df[START_COL].to_numpy()
    ends = df[END_COL].to_numpy()
    _dtype = check_and_return_common_type_2(starts, ends)

    idxs, starts, ends, counts = ruranges.boundary_numpy(  # type: ignore[attr-defined]
        group_ids.astype(np.uint32),
        starts.astype(np.int64),
        ends.astype(np.int64),
    )

    ids = df.take(idxs).loc[:, by]

    result = RangeFrame(
        {START_COL: starts.astype(_dtype), END_COL: ends.astype(_dtype)} | {_by: ids[_by] for _by in by}
    )[col_order]

    return mypy_ensure_rangeframe(result.reset_index(drop=True))


def _outside_bounds(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    df = df.copy()

    _chromsizes = kwargs.get("chromsizes")

    if isinstance(_chromsizes, pd.DataFrame):
        size_df = _chromsizes.df
        if not size_df.Chromosome.is_unique:
            msg = "Chromosomes must be unique in chromsizes."
            raise ValueError(msg)
        chromsizes = dict(*zip(size_df[CHROM_COL], size_df[END_COL], strict=True))
    else:
        if not isinstance(_chromsizes, dict):
            msg = "chromsizes must be a dictionary or a DataFrame"
            raise TypeError(msg)
        chromsizes = _chromsizes

    size = int(chromsizes[df.Chromosome.iloc[0]])
    clip = kwargs.get("clip", False)
    only_right = kwargs.get("only_right", False)

    ends_outright = df.End > size
    starts_outleft = df.Start < 0

    if not clip:  # i.e. remove
        df = df[~ends_outright] if only_right else df[~ends_outright & ~starts_outleft]

    else:
        starts_outright = df.Start >= size

        if only_right:
            df.loc[ends_outright, "End"] = size

            # removing intervals completely out of bounds
            df = df[~starts_outright]

        else:
            ends_outleft = df.End <= 0

            df.loc[ends_outright, "End"] = size
            df.loc[starts_outleft, "Start"] = 0

            # removing intervals completely out of bounds:
            df = df[~starts_outright & ~ends_outleft]

    return df
