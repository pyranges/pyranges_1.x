from typing import TYPE_CHECKING

import pandas as pd

from pyranges.names import CHROM_COL, END_COL, RANGE_COLS, START_COL, STRAND_COL, VALID_STRAND_TYPE
from pyranges.pyranges_helpers import validate_and_convert_strand, mypy_ensure_pyranges

if TYPE_CHECKING:
    from pyranges import PyRanges


def _split(
    df: "PyRanges",
    strand: VALID_STRAND_TYPE = "auto",
    **_,
) -> "PyRanges":
    dtype = df[START_COL].dtype

    starts = df[START_COL]
    ends = df[END_COL]
    points = [starts, ends]

    points = pd.concat(points).sort_values().drop_duplicates()

    _ends = points.shift(-1)

    points = points.iloc[:-1]
    _ends = _ends.iloc[:-1]
    features = pd.concat([points, _ends], axis=1).astype(dtype)
    features.columns = RANGE_COLS

    features.insert(0, CHROM_COL, df[CHROM_COL].iloc[0])
    if validate_and_convert_strand(df, strand):
        features.insert(features.shape[1], STRAND_COL, df[STRAND_COL].iloc[0])

    return mypy_ensure_pyranges(features.reset_index(drop=True))
