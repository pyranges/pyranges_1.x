import numpy as np
import pandas as pd

import pyranges as pr  # noqa: TCH001
from pyranges.core import names
from pyranges.core.example_data import example_data
from pyranges.core.pyranges_helpers import mypy_ensure_pyranges


def random(
    n: int = 1000,
    length: int = 100,
    chromsizes: dict[str, int] | pd.DataFrame | None = None,
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
    0        |    chr11         25516829   25516929   +
    1        |    chr11         132583621  132583721  -
    2        |    chr11         2504795    2504895    +
    3        |    chr11         23816613   23816713   +
    ...      |    ...           ...        ...        ...
    996      |    chr21         30756250   30756350   -
    997      |    chr21         22517078   22517178   +
    998      |    chr21         20605246   20605346   +
    999      |    chr21         21153142   21153242   -
    PyRanges with 1000 rows, 4 columns, and 1 index columns.
    Contains 24 chromosomes and 2 strands.

    """
    rng = np.random.default_rng(seed=seed)

    if chromsizes is None:
        df = example_data.chromsizes
    elif isinstance(chromsizes, dict):
        df = pd.DataFrame({names.CHROM_COL: list(chromsizes.keys()), names.END_COL: list(chromsizes.values())})
    else:
        df = chromsizes

    # Probability of picking each chromosome proportional to its size
    p = df.End / df.End.sum()

    # Determine how many intervals per chromosome
    chosen = rng.choice(df.index, size=n, p=p)
    n_per_chrom = pd.Series(chosen).value_counts(sort=False).to_frame("Count")
    n_per_chrom.insert(1, names.CHROM_COL, pd.Series(df.loc[n_per_chrom.index, names.CHROM_COL].values))

    # Merge chromosome sizes into n_per_chrom for direct access
    n_per_chrom = n_per_chrom.merge(df[[names.CHROM_COL, names.END_COL]], on=names.CHROM_COL, how="left")

    # Extract arrays
    counts_array = n_per_chrom["Count"].to_numpy()
    chroms_array = n_per_chrom[names.CHROM_COL].to_numpy()
    ends_array = n_per_chrom["End"].to_numpy() - length

    # Repeat arrays according to the counts for vectorized generation
    chroms_repeated = np.repeat(chroms_array, counts_array)
    ends_repeated = np.repeat(ends_array, counts_array)

    # Generate random starts in [0, ends_repeated)
    # Using random() gives a uniform [0,1), we scale by ends_repeated
    random_starts = (rng.random(chroms_repeated.size) * ends_repeated).astype(int)

    # Build final DataFrame
    random_df = pd.DataFrame(
        {names.CHROM_COL: chroms_repeated, names.START_COL: random_starts, "End": random_starts + length},
    )

    if strand:
        s = rng.choice(["+", "-"], size=n)
        random_df.insert(3, "Strand", s)

    return mypy_ensure_pyranges(random_df.reset_index(drop=True))
