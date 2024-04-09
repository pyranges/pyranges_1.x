from typing import TYPE_CHECKING

import pandas as pd

from pyranges.core.names import CHROM_COL, END_COL, START_COL, STRAND_COL

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

    res = df.groupby(by, as_index=False).agg(agg_dict)
    return res.reindex(columns=[c for c in col_order if c in res.columns])


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
