import importlib.metadata
import json
import logging
import sys

import numpy as np
import pandas as pd

import pyranges as pr
from pyranges import genomicfeatures, names, statistics
from pyranges.example_data import ExampleData
from pyranges.get_fasta import get_sequence, get_transcript_sequence
from pyranges.methods.concat import concat
from pyranges.multioverlap import count_overlaps
from pyranges.names import END_COL
from pyranges.option_manager import options  # noqa: F401
from pyranges.pyranges_helpers import mypy_ensure_pyranges
from pyranges.pyranges_main import PyRanges  # noqa: F401
from pyranges.range_frame.range_frame import RangeFrame  # noqa: F401
from pyranges.readers import from_string, read_bam, read_bed, read_bigwig, read_gff3, read_gtf  # NOQA: F401

logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


gf = genomicfeatures

__version__ = importlib.metadata.version("pyranges")


TOSTRING_CONSOLE_WIDTH: int | None = None

example_data = ExampleData()
stats = statistics
get_sequence = get_sequence  # noqa: PLW0127
get_transcript_sequence = get_transcript_sequence  # noqa: PLW0127

read_gff = read_gtf

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
        df = pd.DataFrame({names.CHROM_COL: list(chromsizes.keys()), END_COL: list(chromsizes.values())})
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


def version_info() -> None:
    """Print version info for pyranges and its dependencies.

    Used for debugging.
    """
    import importlib.util

    def update_version_info(_version_info: dict[str, str], library: str) -> None:
        version = importlib.import_module(library).__version__ if importlib.util.find_spec(library) else "not installed"

        _version_info[library] = version

    version_info = {
        "pyranges version": pr.__version__,
        "pandas version": pd.__version__,
        "numpy version": np.__version__,
        "python version": ".".join([str(s) for s in sys.version_info]),
    }

    update_version_info(version_info, "ncls")
    update_version_info(version_info, "sorted_nearest")
    update_version_info(version_info, "pyrle")
    update_version_info(version_info, "bamread")
    update_version_info(version_info, "pybigwig")
    update_version_info(version_info, "hypothesis")
    update_version_info(version_info, "pyfaidx")

    LOGGER.info(json.dumps(version_info, indent=4))


__all__ = [
    "from_string",
    "count_overlaps",
    "random",
    "read_gtf",
    "read_bam",
    "read_bed",
    "read_gff3",
    "concat",
    "version_info",
]
