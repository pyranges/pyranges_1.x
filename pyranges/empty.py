import typing
from collections.abc import Iterable

import pandas as pd
from pandas import Series

import pyranges as pr
from pyranges import names

if typing.TYPE_CHECKING:
    from pyranges.pyranges_main import PyRanges


def empty_df(
    with_strand: bool = False,
    columns: Iterable[str] | None = None,
    dtype: Series | None = None,
) -> pd.DataFrame:
    empty = pd.DataFrame(
        columns=list(columns)
        if columns is not None
        else (names.GENOME_LOC_COLS_WITH_STRAND if with_strand else names.GENOME_LOC_COLS),
    )
    return empty.astype(dtype) if dtype is not None else empty


def empty(
    with_strand: bool = False,
    columns: Iterable[str] | None = None,
    dtype: Series | None = None,
) -> "PyRanges":
    """Create an empty PyRanges.

    Parameters
    ----------
    with_strand : bool, default False
        Whether to create a PyRanges with strand information.
    """
    return pr.PyRanges(empty_df(with_strand=with_strand, columns=columns, dtype=dtype))
