from __future__ import print_function

import itertools
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple, Union

import numpy as np
import pandas as pd
import pkg_resources
from natsort import natsorted  # type: ignore

import pyranges as pr
import pyranges.genomicfeatures as gf  # NOQA: F401
from pyranges import data, statistics
from pyranges.get_fasta import get_fasta, get_sequence, get_transcript_sequence
from pyranges.helpers import get_key_from_df, single_value_key
from pyranges.methods.concat import concat
from pyranges.multioverlap import count_overlaps
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
    colnames = ["Chromosome", "Start", "End"]
    if strands is not None:
        if isinstance(strands, str):
            _strands = pd.Series([strands] * len(starts), dtype="category")
        else:
            _strands = pd.Series(strands, dtype="category")

        columns.append(_strands)
        colnames.append("Strand")

    lengths = list(str(len(s)) for s in columns)
    assert len(set(lengths)) == 1, "[{colnames} must be of equal length. But are {columns}".format(
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
    df.columns = pd.Index(colnames)

    return pr.PyRanges(df)


def from_dict(d: Dict[str, Iterable]) -> PyRanges:
    """Create a PyRanges from dict.

    Parameters
    ----------
    d : dict of array-like

        Dict with data.

    Warning
    -------

    On versions of Python prior to 3.6, this function returns a PyRanges with
    the columns in arbitrary order.

    See Also
    --------

    pyranges.from_string : create a PyRanges from a multiline string.

    Examples
    --------

    >>> d = {"Chromosome": [1, 1, 2], "Start": [1, 2, 3], "End": [4, 9, 12], "Strand": ["+", "+", "-"], "ArbitraryValue": ["a", "b", "c"]}
    >>> pr.from_dict(d)
    +--------------+-----------+-----------+--------------+------------------+
    |   Chromosome |     Start |       End | Strand       | ArbitraryValue   |
    |   (category) |   (int64) |   (int64) | (category)   | (object)         |
    |--------------+-----------+-----------+--------------+------------------|
    |            1 |         1 |         4 | +            | a                |
    |            1 |         2 |         9 | +            | b                |
    |            2 |         3 |        12 | -            | c                |
    +--------------+-----------+-----------+--------------+------------------+
    Stranded PyRanges object has 3 rows and 5 columns from 2 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    return PyRanges(pd.DataFrame(d))


def from_dfs(dfs: Union[Dict[str, pd.DataFrame], Dict[Tuple[str, str], pd.DataFrame]]) -> "PyRanges":
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
            raise ValueError(f"All keys must be the same, but df has {_key} and dict had {key}.")

    if not _strand_valid:
        df = pd.concat(empty_removed.values()).reset_index(drop=True)

        groupby_cols = ["Chromosome"]

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

    if strand is None:
        strand = all([gr.stranded for gr in prs])

    if strand is False and any([gr.stranded for gr in prs]):
        prs = [gr.unstrand() for gr in prs]

    grs_per_chromosome = defaultdict(list)
    keys = [gr.dfs.keys() for gr in prs]
    set_keys: Union[Set[str], Set[Tuple[str, str]]] = set(itertools.chain.from_iterable(keys))

    empty_dfs = [pd.DataFrame(columns=gr.columns) for gr in prs]
    for gr, empty in zip(prs, empty_dfs):
        for k in set_keys:
            df = gr.dfs.get(k, empty)  # type: ignore
            grs_per_chromosome[k].append(df)

    if not keys:
        return iter(grs_per_chromosome.values())
    else:
        return iter(natsorted(grs_per_chromosome.items()))


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
        df = pd.DataFrame({"Chromosome": list(chromsizes.keys()), "End": list(chromsizes.values())})
    else:
        df = chromsizes.df

    p = df.End / df.End.sum()

    n_per_chrom = pd.Series(np.random.choice(df.index, size=n, p=p)).value_counts(sort=False).to_frame()
    n_per_chrom.insert(1, "Chromosome", df.loc[n_per_chrom.index].Chromosome)
    n_per_chrom.columns = pd.Index("Count Chromosome".split())

    random_dfs = []
    for _, (count, chrom) in n_per_chrom.iterrows():
        r = np.random.randint(0, df[df.Chromosome == chrom].End - length, size=count)
        _df = pd.DataFrame({"Chromosome": chrom, "Start": r, "End": r + length})
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


def to_bigwig(gr: PyRanges, path: Path, chromosome_sizes=Optional[Chromsizes]):
    """Write df to bigwig.

    Must contain the columns Chromosome, Start, End and Score. All others are ignored.

    Parameters
    ----------
    gr: PyRanges
        Intervals to write.

    path : Path

        Where to write bigwig.

    chromosome_sizes : PyRanges or dict

        If dict: map of chromosome names to chromosome length.

    Examples
    --------

    Extended example with how to prepare your data for writing bigwigs:

    >>> d =  {'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [1, 4, 6],
    ...       'End': [7, 8, 10], 'Strand': ['+', '-', '-'],
    ...       'Value': [10, 20, 30]}
    >>> import pyranges as pr
    >>> gr = pr.from_dict(d)
    >>> hg19 = pr.data.chromsizes()
    >>> print(hg19)
    +--------------+-----------+-----------+
    | Chromosome   | Start     | End       |
    | (category)   | (int64)   | (int64)   |
    |--------------+-----------+-----------|
    | chr1         | 0         | 249250621 |
    | chr2         | 0         | 243199373 |
    | chr3         | 0         | 198022430 |
    | chr4         | 0         | 191154276 |
    | ...          | ...       | ...       |
    | chr22        | 0         | 51304566  |
    | chrM         | 0         | 16571     |
    | chrX         | 0         | 155270560 |
    | chrY         | 0         | 59373566  |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 25 rows and 3 columns from 25 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.

    Overlapping intervals are invalid in bigwigs:

    >>> to_bigwig(gr, "outpath.bw", hg19)
    Traceback (most recent call last):
    ...
    AssertionError: Can only write one strand at a time. Use an unstranded PyRanges or subset on strand first.

    >>> to_bigwig(gr["-"], "outpath.bw", hg19)
    Traceback (most recent call last):
    ...
    AssertionError: Intervals must not overlap.

    >>> gr
    +--------------+-----------+-----------+--------------+-----------+
    | Chromosome   |     Start |       End | Strand       |     Value |
    | (category)   |   (int64) |   (int64) | (category)   |   (int64) |
    |--------------+-----------+-----------+--------------+-----------|
    | chr1         |         1 |         7 | +            |        10 |
    | chr1         |         4 |         8 | -            |        20 |
    | chr1         |         6 |        10 | -            |        30 |
    +--------------+-----------+-----------+--------------+-----------+
    Stranded PyRanges object has 3 rows and 5 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    >>> value = gr.to_rle(rpm=False, value_col="Value")
    >>> value
    chr1 +
    --
    +--------+-----+------+
    | Runs   | 1   | 6    |
    |--------+-----+------|
    | Values | 0.0 | 10.0 |
    +--------+-----+------+
    Rle of length 7 containing 2 elements (avg. length 3.5)
    <BLANKLINE>
    chr1 -
    --
    +--------+-----+------+------+------+
    | Runs   | 4   | 2    | 2    | 2    |
    |--------+-----+------+------+------|
    | Values | 0.0 | 20.0 | 50.0 | 30.0 |
    +--------+-----+------+------+------+
    Rle of length 10 containing 4 elements (avg. length 2.5)
    RleDict object with 2 chromosomes/strand pairs.

    >>> raw = gr.to_rle(rpm=False)
    >>> raw
    chr1 +
    --
    +--------+-----+-----+
    | Runs   | 1   | 6   |
    |--------+-----+-----|
    | Values | 0.0 | 1.0 |
    +--------+-----+-----+
    Rle of length 7 containing 2 elements (avg. length 3.5)
    <BLANKLINE>
    chr1 -
    --
    +--------+-----+-----+-----+-----+
    | Runs   | 4   | 2   | 2   | 2   |
    |--------+-----+-----+-----+-----|
    | Values | 0.0 | 1.0 | 2.0 | 1.0 |
    +--------+-----+-----+-----+-----+
    Rle of length 10 containing 4 elements (avg. length 2.5)
    RleDict object with 2 chromosomes/strand pairs.

    >>> result = (value / raw).apply_values(np.log10)
    >>> result
    chr1 +
    --
    +--------+-----+-----+
    | Runs   | 1   | 6   |
    |--------+-----+-----|
    | Values | nan | 1.0 |
    +--------+-----+-----+
    Rle of length 7 containing 2 elements (avg. length 3.5)
    <BLANKLINE>
    chr1 -
    --
    +--------+-----+--------------------+--------------------+--------------------+
    | Runs   | 4   | 2                  | 2                  | 2                  |
    |--------+-----+--------------------+--------------------+--------------------|
    | Values | nan | 1.3010300397872925 | 1.3979400396347046 | 1.4771212339401245 |
    +--------+-----+--------------------+--------------------+--------------------+
    Rle of length 10 containing 4 elements (avg. length 2.5)
    RleDict object with 2 chromosomes/strand pairs.

    >>> out = result.numbers_only().to_ranges()
    >>> out
    +--------------+-----------+-----------+-------------+--------------+
    | Chromosome   |     Start |       End |       Score | Strand       |
    | (category)   |   (int64) |   (int64) |   (float64) | (category)   |
    |--------------+-----------+-----------+-------------+--------------|
    | chr1         |         1 |         7 |     1       | +            |
    | chr1         |         4 |         6 |     1.30103 | -            |
    | chr1         |         6 |         8 |     1.39794 | -            |
    | chr1         |         8 |        10 |     1.47712 | -            |
    +--------------+-----------+-----------+-------------+--------------+
    Stranded PyRanges object has 4 rows and 5 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    >>> to_bigwig(out["-"], "deleteme_reverse.bw", hg19)
    >>> to_bigwig(out["+"], "deleteme_forward.bw", hg19)
    """

    try:
        import pyBigWig  # type: ignore
    except ModuleNotFoundError:
        print(
            "pybigwig must be installed to create bigwigs. Use `conda install -c bioconda pybigwig` or `pip install pybigwig` to install it."
        )
        import sys

        sys.exit(1)

    assert (
        len(gr.strands) <= 1
    ), "Can only write one strand at a time. Use an unstranded PyRanges or subset on strand first."
    lengths = gr.lengths()
    assert isinstance(lengths, pd.Series)
    assert np.sum(lengths) == gr.merge().length, "Intervals must not overlap."

    df = gr.df

    unique_chromosomes = list(df.Chromosome.drop_duplicates())

    if not isinstance(chromosome_sizes, dict):
        size_df = chromosome_sizes.df
        chromosome_sizes = {k: v for k, v in zip(size_df.Chromosome, size_df.End)}

    header = [(c, int(chromosome_sizes[c])) for c in unique_chromosomes]

    bw = pyBigWig.open(path, "w")
    bw.addHeader(header)

    chromosomes = df.Chromosome.tolist()
    starts = df.Start.tolist()
    ends = df.End.tolist()
    values = df.Score.tolist()

    bw.addEntries(chromosomes, starts, ends=ends, values=values)


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
    "to_bigwig",
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


def _test():
    import doctest

    doctest.testmod()


if __name__ == "__main__":
    _test()
