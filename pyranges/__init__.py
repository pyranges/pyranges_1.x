from __future__ import print_function

from pyranges.range_frame import RangeFrame

import itertools
import sys
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple, Union

import numpy as np
import pandas as pd
import pkg_resources
from natsort import natsorted  # type: ignore
from pandas import Series

import pyranges as pr
import pyranges.genomicfeatures as gf  # NOQA: F401
from pyranges import data, statistics
from pyranges.get_fasta import get_fasta, get_sequence, get_transcript_sequence
from pyranges.helpers import get_key_from_df, single_value_key
from pyranges.methods.concat import concat
from pyranges.multioverlap import count_overlaps
from pyranges.names import (
    GENOME_LOC_COLS_WITH_STRAND,
    GENOME_LOC_COLS,
    STRAND_COL,
    CHROM_COL,
    START_COL,
    END_COL,
)
from pyranges.pyranges_main import PyRanges
from pyranges.readers import read_bam, read_bed, read_bigwig, read_gff3, read_gtf  # NOQA: F401

__version__ = pkg_resources.get_distribution("pyranges").version

get_example_path = data.get_example_path
stats = statistics
get_fasta = get_fasta
get_sequence = get_sequence
get_transcript_sequence = get_transcript_sequence

read_gff = read_gtf

Chromsizes = Union[Dict[str, int], Dict[Tuple[str, str], int]]


def concat(grs: Iterable[PyRanges], *args, **kwargs) -> "PyRanges":
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
    >>> gr1 = pr.data.f2()
    >>> gr2 = pr.data.f1()
    >>> pr.concat([gr1, gr2])
    Chromosome      Start      End  Name         Score  Strand
    object          int64    int64  object       int64  object
    ------------  -------  -------  ---------  -------  --------
    chr1                1        2  a                0  +
    chr1                6        7  b                0  -
    chr1                3        6  interval1        0  +
    chr1                5        7  interval2        0  -
    chr1                8        9  interval3        0  +
    PyRanges with 5 rows and 6 columns.
    Contains 1 chromosomes and 2 strands.

    >>> pr.concat([gr1, gr2.remove_strand()])
    Chromosome      Start      End  Name         Score  Strand
    object          int64    int64  object       int64  object
    ------------  -------  -------  ---------  -------  --------
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
    columns: Optional[Iterable[str]] = None,
    dtype: Optional[Series] = None,
) -> pd.DataFrame:
    empty = pd.DataFrame(
        columns=columns
        if list(columns)
        else (GENOME_LOC_COLS_WITH_STRAND if with_strand else GENOME_LOC_COLS),
    )
    return empty.astype(dtype) if dtype is not None else empty


def empty(
    with_strand: bool = False,
    columns: Optional[Iterable[str]] = None,
    dtype: Optional[Series] = None,
) -> "PyRanges":
    """Create an empty PyRanges.

    Parameters
    ----------
    with_strand : bool, default False
        Whether to create a PyRanges with strand information.
    """
    return pr.PyRanges(empty_df(with_strand=with_strand, columns=columns, dtype=dtype))


def from_args(
    chromosomes: Union[Sequence[str], Sequence[int]],
    starts: Sequence[int],
    ends: Sequence[int],
    strands: Optional[Union[str, Sequence[str]]] = None,
) -> "PyRanges":
    if isinstance(chromosomes, str) or isinstance(chromosomes, int):
        _chromosomes = pd.Series([chromosomes] * len(starts), dtype="category")
    else:
        _chromosomes = pd.Series(chromosomes, dtype="category")

    columns: List[pd.Series] = [_chromosomes, pd.Series(starts), pd.Series(ends)]
    colnames = GENOME_LOC_COLS[:]
    if strands is not None:
        if isinstance(strands, str):
            _strands = pd.Series([strands] * len(starts), dtype="category")
        else:
            _strands = pd.Series(strands, dtype="category")

        columns.append(_strands)
        colnames.append(STRAND_COL)

    lengths = list(str(len(s)) for s in columns)
    assert (
        len(set(lengths)) == 1
    ), "[{colnames} must be of equal length. But are {columns}".format(
        colnames=", ".join(colnames), columns=", ".join(lengths)
    )

    idx = range(len(starts))
    series_to_concat = []
    for s in columns:
        if isinstance(s, pd.Series):
            s = pd.Series(s.values, index=idx)
        else:
            s = pd.Series(s, index=idx)

        series_to_concat.append(s)

    df = pd.concat(series_to_concat, axis=1)
    df.columns = colnames

    return pr.PyRanges(df)


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


def from_string(s: str) -> PyRanges:
    """Create a PyRanges from multiline string.

    Parameters
    ----------
    s : str

        String with data.

    See Also
    --------

    pyranges.from_dict : create a PyRanges from a dictionary.

    Examples
    --------

    >>> s = '''Chromosome      Start        End Strand
    ... chr1  246719402  246719502      +
    ... chr5   15400908   15401008      +
    ... chr9   68366534   68366634      +
    ... chr14   79220091   79220191      +
    ... chr14  103456471  103456571      -'''

    >>> pr.from_string(s)
    +--------------+-----------+-----------+--------------+
    | Chromosome   |     Start |       End | Strand       |
    | (category)   |   (int64) |   (int64) | (category)   |
    |--------------+-----------+-----------+--------------|
    | chr1         | 246719402 | 246719502 | +            |
    | chr5         |  15400908 |  15401008 | +            |
    | chr9         |  68366534 |  68366634 | +            |
    | chr14        |  79220091 |  79220191 | +            |
    | chr14        | 103456471 | 103456571 | -            |
    +--------------+-----------+-----------+--------------+
    Stranded PyRanges object has 5 rows and 4 columns from 4 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    """

    from io import StringIO

    df = pd.read_csv(StringIO(s), sep=r"\s+", index_col=None)

    return PyRanges(df)


def itergrs(prs: Iterable[PyRanges], strand=None, keys=False):
    r"""Iterate over multiple PyRanges at once.

    Parameters
    ----------
    prs : list of PyRanges

        PyRanges to iterate over.

    strand : bool, default None, i.e. auto

        Whether to iterate over strands. If True, all PyRanges must be stranded.

    keys : bool, default False

        Return tuple with key and value from iterator.

    Examples
    --------

    >>> d1 = {"Chromosome": [1, 1, 2], "Start": [1, 2, 3], "End": [4, 9, 12], "Strand": ["+", "+", "-"]}
    >>> d2 = {"Chromosome": [2, 3, 3], "Start": [5, 9, 21], "End": [81, 42, 25], "Strand": ["-", "+", "-"]}
    >>> gr1, gr2 = pr.from_dict(d1), pr.from_dict(d2)
    >>> gr1
    +--------------+-----------+-----------+--------------+
    |   Chromosome |     Start |       End | Strand       |
    |   (category) |   (int64) |   (int64) | (category)   |
    |--------------+-----------+-----------+--------------|
    |            1 |         1 |         4 | +            |
    |            1 |         2 |         9 | +            |
    |            2 |         3 |        12 | -            |
    +--------------+-----------+-----------+--------------+
    Stranded PyRanges object has 3 rows and 4 columns from 2 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    >>> gr2
    +--------------+-----------+-----------+--------------+
    |   Chromosome |     Start |       End | Strand       |
    |   (category) |   (int64) |   (int64) | (category)   |
    |--------------+-----------+-----------+--------------|
    |            2 |         5 |        81 | -            |
    |            3 |         9 |        42 | +            |
    |            3 |        21 |        25 | -            |
    +--------------+-----------+-----------+--------------+
    Stranded PyRanges object has 3 rows and 4 columns from 2 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    >>> ranges = [gr1, gr2]

    >>> for key, dfs in pr.itergrs(ranges, keys=True):
    ...     print("-----------\n" + str(key) + "\n-----------")
    ...     for df in dfs:
    ...         print(df)
    -----------
    ('1', '+')
    -----------
      Chromosome  Start  End Strand
    0          1      1    4      +
    1          1      2    9      +
    Empty DataFrame
    Columns: [Chromosome, Start, End, Strand]
    Index: []
    -----------
    ('2', '-')
    -----------
      Chromosome  Start  End Strand
    2          2      3   12      -
      Chromosome  Start  End Strand
    0          2      5   81      -
    -----------
    ('3', '+')
    -----------
    Empty DataFrame
    Columns: [Chromosome, Start, End, Strand]
    Index: []
      Chromosome  Start  End Strand
    1          3      9   42      +
    -----------
    ('3', '-')
    -----------
    Empty DataFrame
    Columns: [Chromosome, Start, End, Strand]
    Index: []
      Chromosome  Start  End Strand
    2          3     21   25      -
    """

    # Determine if all prs are stranded only if strand is None, otherwise use the provided value.
    strand = all(gr.strand_values_valid for gr in prs) if strand is None else strand

    # If strand is False and any PyRanges are stranded, remove strand information.
    prs = [gr.remove_strand() for gr in prs if gr.strand_values_valid] if not strand else prs

    empty_dfs = [pd.DataFrame(columns=gr.columns) for gr in prs]
    for gr, empty in zip(prs, empty_dfs):
        for k in set_keys:
            df = gr.dfs.get(k, empty)  # type: ignore
            grs_per_chromosome[k].append(df)

    # Iterate through each key and PyRanges object, collecting DataFrames.
    for key in set_keys:
        yield key, [gr.dfs.get(key, pd.DataFrame(columns=gr.columns)) for gr in prs]


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

    # >>> pr.random()
    # +--------------+-----------+-----------+--------------+
    # | Chromosome   | Start     | End       | Strand       |
    # | (category)   | (int64)   | (int64)   | (category)   |
    # |--------------+-----------+-----------+--------------|
    # | chr1         | 216128004 | 216128104 | +            |
    # | chr1         | 114387955 | 114388055 | +            |
    # | chr1         | 67597551  | 67597651  | +            |
    # | chr1         | 26306616  | 26306716  | +            |
    # | ...          | ...       | ...       | ...          |
    # | chrY         | 20811459  | 20811559  | -            |
    # | chrY         | 12221362  | 12221462  | -            |
    # | chrY         | 8578041   | 8578141   | -            |
    # | chrY         | 43259695  | 43259795  | -            |
    # +--------------+-----------+-----------+--------------+
    # Stranded PyRanges object has 1,000 rows and 4 columns from 24 chromosomes.
    # For printing, the PyRanges was sorted on Chromosome and Strand.

    To have random interval lengths:

    # >>> gr = pr.random(length=1)
    # >>> gr.End += np.random.randint(int(1e5), size=len(gr))
    # >>> gr.Length = gr.lengths()
    # >>> gr
    # +--------------+-----------+-----------+--------------+-----------+
    # | Chromosome   | Start     | End       | Strand       | Length    |
    # | (category)   | (int64)   | (int64)   | (category)   | (int64)   |
    # |--------------+-----------+-----------+--------------+-----------|
    # | chr1         | 203654331 | 203695380 | +            | 41049     |
    # | chr1         | 46918271  | 46978908  | +            | 60637     |
    # | chr1         | 97355021  | 97391587  | +            | 36566     |
    # | chr1         | 57284999  | 57323542  | +            | 38543     |
    # | ...          | ...       | ...       | ...          | ...       |
    # | chrY         | 31665821  | 31692660  | -            | 26839     |
    # | chrY         | 20236607  | 20253473  | -            | 16866     |
    # | chrY         | 33255377  | 33315933  | -            | 60556     |
    # | chrY         | 31182964  | 31205467  | -            | 22503     |
    # +--------------+-----------+-----------+--------------+-----------+
    # Stranded PyRanges object has 1,000 rows and 5 columns from 24 chromosomes.
    # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    if chromsizes is None:
        df = data.chromsizes().df
    elif isinstance(chromsizes, dict):
        df = pd.DataFrame(
            {CHROM_COL: list(chromsizes.keys()), "End": list(chromsizes.values())}
        )
    else:
        df = chromsizes.df

    p = df.End / df.End.sum()

    n_per_chrom = (
        pd.Series(np.random.choice(df.index, size=n, p=p))
        .value_counts(sort=False)
        .to_frame()
    )
    n_per_chrom.insert(1, CHROM_COL, df.loc[n_per_chrom.index].Chromosome)
    n_per_chrom.columns = pd.Index("Count Chromosome".split())

    random_dfs = []
    for _, (count, chrom) in n_per_chrom.iterrows():
        r = np.random.randint(0, df[df.Chromosome == chrom].End - length, size=count)
        _df = pd.DataFrame({CHROM_COL: chrom, START_COL: r, "End": r + length})
        random_dfs.append(_df)

    random_df = pd.concat(random_dfs)
    if strand:
        s = np.random.choice("+ -".split(), size=n)
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
    update_version_info(version_info, "ray")
    update_version_info(version_info, "bamread")
    # update_version_info(version_info, "bwread") no version string yet!
    update_version_info(version_info, "pyranges_db")
    update_version_info(version_info, "pybigwig")
    update_version_info(version_info, "hypothesis")

    print(version_info)


__all__ = [
    "from_string",
    "from_dict",
    "count_overlaps",
    "random",
    "itergrs",
    "read_gtf",
    "read_bam",
    "read_bed",
    "read_gff3",
    "concat",
    "PyRanges",
    "version_info",
]
