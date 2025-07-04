from typing import TYPE_CHECKING

import pandas as pd

from pyranges.core.names import CHROM_COL, END_COL, START_COL
from pyranges.core.pyranges_helpers import ensure_pyranges

if TYPE_CHECKING:
    from pyranges import PyRanges


def tile_genome(
    chromsizes: "PyRanges | pd.DataFrame | dict[str | int, int]",
    tile_size: int,
    *,
    full_last_tile: bool = False,
) -> "PyRanges":
    """Create a tiled genome.

    Split the genome into adjacent non-overlapping tiles of a given size.

    Parameters
    ----------
    chromsizes : dict or PyRanges
        Dict or PyRanges describing the lengths of the chromosomes.

    tile_size : int
        Length of the tiles.

    full_last_tile : bool, default False
        Do not truncate the last tile to the end of the chromosome. Use to ensure size is consistent for all tiles.

    See Also
    --------
    pyranges.PyRanges.tile_ranges : split intervals into adjacent non-overlapping tiles.

    Examples
    --------
    >>> import pyranges as pr
    >>> chromsizes = pr.example_data.chromsizes
    >>> chromsizes
    index    |    Chromosome    Start    End
    int64    |    category      int64    int64
    -------  ---  ------------  -------  ---------
    0        |    chr1          0        249250621
    1        |    chr2          0        243199373
    2        |    chr3          0        198022430
    3        |    chr4          0        191154276
    ...      |    ...           ...      ...
    21       |    chr19         0        59128983
    22       |    chr22         0        51304566
    23       |    chr21         0        48129895
    24       |    chrM          0        16571
    PyRanges with 25 rows, 3 columns, and 1 index columns.
    Contains 25 chromosomes.

    >>> pr.tile_genome(chromsizes,int(1e6))
    index    |    Chromosome    Start     End
    int64    |    category      int64     int64
    -------  ---  ------------  --------  --------
    0        |    chr1          0         1000000
    1        |    chr1          1000000   2000000
    2        |    chr1          2000000   3000000
    3        |    chr1          3000000   4000000
    ...      |    ...           ...       ...
    3110     |    chrY          56000000  57000000
    3111     |    chrY          57000000  58000000
    3112     |    chrY          58000000  59000000
    3113     |    chrY          59000000  59373566
    PyRanges with 3114 rows, 3 columns, and 1 index columns.
    Contains 25 chromosomes.

    >>> pr.tile_genome(chromsizes,int(1e6), full_last_tile=True)
    index    |    Chromosome    Start     End
    int64    |    category      int64     int64
    -------  ---  ------------  --------  --------
    0        |    chr1          0         1000000
    1        |    chr1          1000000   2000000
    2        |    chr1          2000000   3000000
    3        |    chr1          3000000   4000000
    ...      |    ...           ...       ...
    3110     |    chrY          56000000  57000000
    3111     |    chrY          57000000  58000000
    3112     |    chrY          58000000  59000000
    3113     |    chrY          59000000  60000000
    PyRanges with 3114 rows, 3 columns, and 1 index columns.
    Contains 25 chromosomes.

    """
    if isinstance(chromsizes, dict):
        chromsize_dict = chromsizes
        chromosomes, ends = list(chromsizes.keys()), list(chromsizes.values())
        df = pd.DataFrame({CHROM_COL: chromosomes, START_COL: 0, END_COL: ends})
        chromsizes = pd.DataFrame(df)
    else:
        chromsize_dict = dict(zip(chromsizes[CHROM_COL], chromsizes[END_COL], strict=True))

    gr = ensure_pyranges(chromsizes).sort_ranges().tile_ranges(tile_size)

    return ensure_pyranges(
        gr.reset_index(drop=True)
        if full_last_tile
        else gr.groupby(CHROM_COL).apply(_last_tile, sizes=chromsize_dict).reset_index(drop=True),
    )


def _last_tile(df: pd.DataFrame, sizes: dict[str, int]) -> pd.DataFrame:
    size = sizes[df.Chromosome.iloc[0]]
    df.iloc[-1, [*df.columns].index(END_COL)] = size
    return df
