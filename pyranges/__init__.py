from __future__ import print_function

TOSTRING_CONSOLE_WIDTH = None

import sys
from typing import Dict, Iterable, Optional, Tuple, Union, TYPE_CHECKING

import numpy as np
import pandas as pd
import pkg_resources
from pandas import Series
from pyranges.range_frame import RangeFrame
import pyranges as pr
import pyranges.genomicfeatures as gf
from pyranges.pyranges_main import PyRanges
from pyranges.example_data import ExampleData
from pyranges import data, statistics
from pyranges.get_fasta import get_fasta, get_sequence, get_transcript_sequence
from pyranges.helpers import get_key_from_df, single_value_key
from pyranges.methods.concat import concat
from pyranges.multioverlap import count_overlaps
from pyranges.names import (
    GENOME_LOC_COLS_WITH_STRAND,
    GENOME_LOC_COLS,
    CHROM_COL,
    START_COL,
)
if TYPE_CHECKING:
    from pyranges import PyRanges
from pyranges.readers import read_bam, read_bed, read_bigwig, read_gff3, read_gtf  # NOQA: F401

__version__ = pkg_resources.get_distribution("pyranges").version

data = ExampleData()
stats = statistics
get_fasta = get_fasta
get_sequence = get_sequence
get_transcript_sequence = get_transcript_sequence

read_gff = read_gtf

Chromsizes = Union[Dict[str, int], Dict[Tuple[str, str], int]]


def concat(grs: Iterable["PyRanges"], *args, **kwargs) -> "PyRanges":
    """
    Concatenate PyRanges.

    Parameters
    ----------
    grs

    Returns
    -------
    pyranges.PyRanges

    Examples:
    ---------
    >>> gr1 = pr.data.f2
    >>> gr2 = pr.data.f1
    >>> pr.concat([gr1, gr2])
    Chromosome      Start      End  Name         Score  Strand
    category        int64    int64  object       int64  category
    ------------  -------  -------  ---------  -------  ----------
    chr1                1        2  a                0  +
    chr1                6        7  b                0  -
    chr1                3        6  interval1        0  +
    chr1                5        7  interval2        0  -
    chr1                8        9  interval3        0  +
    PyRanges with 5 rows and 6 columns.
    Contains 1 chromosomes and 2 strands.

    >>> pr.concat([gr1, gr2.remove_strand()])
    Chromosome      Start      End  Name         Score  Strand
    category        int64    int64  object       int64  category
    ------------  -------  -------  ---------  -------  ----------
    chr1                1        2  a                0  +
    chr1                6        7  b                0  -
    chr1                3        6  interval1        0  nan
    chr1                5        7  interval2        0  nan
    chr1                8        9  interval3        0  nan
    PyRanges with 5 rows and 6 columns.
    Contains 1 chromosomes and 2 strands (including non-genomic strands: nan).
    """
    return PyRanges(pd.concat(grs, *args, **kwargs))


def empty_df(
    with_strand: bool = False,
    columns: Iterable[str] | None = None,
    dtype: Series | None = None,
) -> pd.DataFrame:
    empty = pd.DataFrame(
        columns=list(columns)
        if columns is not None
        else (GENOME_LOC_COLS_WITH_STRAND if with_strand else GENOME_LOC_COLS),
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


def from_dfs(
    dfs: Union[Dict[str, pd.DataFrame], Dict[Tuple[str, str], pd.DataFrame]]
) -> "PyRanges":
    df: pd.DataFrame
    empty_removed = {k: v.copy() for k, v in dfs.items() if not v.empty}

    _strand_valid = True
    for key, df in empty_removed.items():
        _key = get_key_from_df(df)
        if not single_value_key(df):
            raise ValueError("All Chromosome/Strand vals in a df must be the same.")
        _key_same = _key == key

        if isinstance(_key, tuple):
            _strand_valid = _strand_valid and (_key[1] in ["+", "-"])

        if _strand_valid and not _key_same:
            raise ValueError(
                f"All keys must be the same, but df has {_key} and dict had {key}."
            )

    if not _strand_valid:
        df = pd.concat(empty_removed.values()).reset_index(drop=True)

        groupby_cols = [CHROM_COL]

        empty_removed = {k[0]: v for k, v in df.groupby(groupby_cols)}  # type: ignore

    gr = PyRanges()
    gr.__dict__["dfs"] = empty_removed

    return gr  # type: ignore


def from_string(s: str) -> "PyRanges":
    """Create a PyRanges from multiline string.

    Parameters
    ----------
    s : str

        String with data.

    Examples
    --------

    >>> s = '''Chromosome      Start        End Strand
    ... chr1  246719402  246719502      +
    ... chr5   15400908   15401008      +
    ... chr9   68366534   68366634      +
    ... chr14   79220091   79220191      +
    ... chr14  103456471  103456571      -'''

    >>> pr.from_string(s)
    Chromosome        Start        End  Strand
    object            int64      int64  object
    ------------  ---------  ---------  --------
    chr1          246719402  246719502  +
    chr5           15400908   15401008  +
    chr9           68366534   68366634  +
    chr14          79220091   79220191  +
    chr14         103456471  103456571  -
    PyRanges with 5 rows and 4 columns.
    Contains 4 chromosomes and 2 strands.
    """

    from io import StringIO

    df = pd.read_csv(StringIO(s), sep=r"\s+", index_col=None)

    return PyRanges(df)


def random(
    n: int = 1000,
    length: int = 100,
    chromsizes: Optional[Chromsizes] = None,
    strand: bool = True,
    seed: Optional[int] = None,
):
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
    Chromosome    Start      End        Strand
    object        int64      int64      object
    ------------  ---------  ---------  --------
    chr4          130788360  130788460  +
    chr4          36129012   36129112   +
    chr4          69733790   69733890   -
    chr4          187723767  187723867  -
    ...           ...        ...        ...
    chr21         13544178   13544278   -
    chr21         33556472   33556572   +
    chr21         31438477   31438577   +
    chr21         38433522   38433622   -
    PyRanges with 1000 rows and 4 columns.
    Contains 24 chromosomes and 2 strands.
    """
    rng = np.random.default_rng(seed=seed)

    if chromsizes is None:
        df = data.chromsizes
    elif isinstance(chromsizes, dict):
        df = pd.DataFrame(
            {CHROM_COL: list(chromsizes.keys()), "End": list(chromsizes.values())}
        )
    else:
        df = chromsizes

    p = df.End / df.End.sum()

    n_per_chrom = (
        pd.Series(rng.choice(df.index, size=n, p=p)).value_counts(sort=False).to_frame()
    )
    n_per_chrom.insert(1, CHROM_COL, df.loc[n_per_chrom.index].Chromosome)
    n_per_chrom.columns = pd.Index("Count Chromosome".split())

    random_dfs = []
    for _, (count, chrom) in n_per_chrom.iterrows():
        r = rng.integers(0, df[df.Chromosome == chrom].End - length, size=count)
        _df = pd.DataFrame({CHROM_COL: chrom, START_COL: r, "End": r + length})
        random_dfs.append(_df)

    random_df = pd.concat(random_dfs)
    if strand:
        s = rng.choice("+ -".split(), size=n)
        random_df.insert(3, "Strand", s)

    return PyRanges(random_df)


"""Namespace for statistcal functions.

See Also
--------
pyranges.statistics : statistcal methods for genomics."""


def version_info() -> None:
    import importlib

    def update_version_info(_version_info, library) -> None:
        if importlib.util.find_spec(library):  # type: ignore
            version = importlib.import_module(library).__version__
        else:
            version = "not installed"

        _version_info[library] = version

    version_info = {
        "pyranges version": pr.__version__,
        "pandas version": pd.__version__,
        "numpy version": np.__version__,
        "python version": sys.version_info,
    }

    update_version_info(version_info, "ncls")
    update_version_info(version_info, "sorted_nearest")
    update_version_info(version_info, "pyrle")
    update_version_info(version_info, "bamread")
    update_version_info(version_info, "pyranges_db")
    update_version_info(version_info, "pybigwig")
    update_version_info(version_info, "hypothesis")

    print(version_info)


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
