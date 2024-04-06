import numpy as np
import pandas as pd

import pyranges as pr
from pyranges import names
from pyranges.example_data_manager import example_data
from pyranges.pyranges_helpers import mypy_ensure_pyranges

Chromsizes = dict[str, int] | dict[tuple[str, str], int]


def random(
    n: int = 1000,
    length: int = 100,
    chromsizes: Chromsizes | None = None,
    seed: int | None = None,
    *,
    strand: bool = True,
) -> "pr.PyRanges":
    """Return PyRanges with random intervals.

    Parameters
    ----------
    n : int, default 1000
        Number of intervals.

    length : int, default 100
        Length of intervals.

    chromsizes : dict or DataFrame, default None, i.e. use "hg19"
        Draw intervals from within these bounds.

    strand : bool, default True
        Data should have strand.

    seed : int, default None
        Seed for random number generator.

    Examples
    --------
    >>> pr.random(seed=12345)
    index    |    Chromosome    Start      End        Strand
    int64    |    object        int64      int64      object
    -------  ---  ------------  ---------  ---------  --------
    0        |    chr4          130788360  130788460  +
    1        |    chr4          36129012   36129112   +
    2        |    chr4          69733790   69733890   -
    3        |    chr4          187723767  187723867  -
    ...      |    ...           ...        ...        ...
    996      |    chr21         13544178   13544278   -
    997      |    chr21         33556472   33556572   +
    998      |    chr21         31438477   31438577   +
    999      |    chr21         38433522   38433622   -
    PyRanges with 1000 rows, 4 columns, and 1 index columns.
    Contains 24 chromosomes and 2 strands.

    """
    rng = np.random.default_rng(seed=seed)

    df: pd.DataFrame
    if chromsizes is None:
        df = example_data.chromsizes
    elif isinstance(chromsizes, dict):
        df = pd.DataFrame({names.CHROM_COL: list(chromsizes.keys()), names.END_COL: list(chromsizes.values())})
    else:
        df = chromsizes

    p = df.End / df.End.sum()

    n_per_chrom = pd.Series(rng.choice(df.index, size=n, p=p)).value_counts(sort=False).to_frame()
    n_per_chrom.insert(1, names.CHROM_COL, df.loc[n_per_chrom.index].Chromosome)
    n_per_chrom.columns = pd.Index("Count Chromosome".split())

    random_dfs = []
    for _, (count, chrom) in n_per_chrom.iterrows():
        r = rng.integers(0, df[df.Chromosome == chrom].End - length, size=count)
        _df = pd.DataFrame({names.CHROM_COL: chrom, names.START_COL: r, "End": r + length})
        random_dfs.append(_df)

    random_df = pd.concat(random_dfs)
    if strand:
        s = rng.choice("+ -".split(), size=n)
        random_df.insert(3, "Strand", s)

    return mypy_ensure_pyranges(random_df.reset_index(drop=True))
