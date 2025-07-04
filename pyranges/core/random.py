from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyranges.core import names
from pyranges.core.example_data import example_data
from pyranges.core.pyranges_helpers import ensure_pyranges

if TYPE_CHECKING:
    import pyranges as pr


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
    >>> import pyranges as pr
    >>> pr.random(seed=12345)
    index    |    Chromosome    Start      End        Strand
    int64    |    object        int64      int64      object
    -------  ---  ------------  ---------  ---------  --------
    0        |    chr4          36129012   36129112   +
    1        |    chr5          177668498  177668598  -
    2        |    chr15         1902279    1902379    +
    3        |    chr11         23816613   23816713   +
    ...      |    ...           ...        ...        ...
    996      |    chr2          155410960  155411060  -
    997      |    chr6          80054552   80054652   +
    998      |    chrX          66474125   66474225   +
    999      |    chr7          69941721   69941821   -
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

    # Pick chromosomes (by picking rows from df)
    chosen = rng.choice(df.index, size=n, p=p)

    # Extract the chosen chromosomes and their lengths
    chrom_chosen = df.loc[chosen, names.CHROM_COL].to_numpy()
    end_chosen = df.loc[chosen, names.END_COL].to_numpy() - length

    # Generate random starts, ensuring they don't exceed (End - length)
    random_starts = (rng.random(n) * end_chosen).astype(int)

    # Construct final DataFrame
    random_df = pd.DataFrame(
        {names.CHROM_COL: chrom_chosen, names.START_COL: random_starts, "End": random_starts + length}
    )

    if strand:
        s = rng.choice(["+", "-"], size=n)
        random_df["Strand"] = s

    return ensure_pyranges(random_df.reset_index(drop=True))
