from typing import TYPE_CHECKING, Optional

from pyranges.methods import concat
from pyranges.names import GENOME_LOC_COLS_WITH_STRAND, VALID_OVERLAP_TYPE, VALID_STRAND_BEHAVIOR_TYPE

if TYPE_CHECKING:
    from pyranges.pyranges_main import PyRanges


def count_overlaps(
    grs: dict[str, "PyRanges"],
    features: Optional["PyRanges"] = None,
    strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
    by: list[str] | None = None,
    how: VALID_OVERLAP_TYPE = "all",
) -> "PyRanges":
    """Count overlaps in multiple pyranges.

    Parameters
    ----------
    grs : dict of PyRanges

        The PyRanges to use as queries.

    features : PyRanges, default None

        The PyRanges to use as subject in the query. If None, the PyRanges themselves are used as a query.

    strand_behavior : {None, "same", "opposite", False}, default None, i.e. auto

        Whether to compare PyRanges on the same strand, the opposite or ignore strand
        information. The default, None, means use "same" if both PyRanges are stranded,
        otherwise ignore the strand information.

     how : {None, "all", "containment"}, default None, i.e. all

        What intervals to report. By default reports all overlapping intervals. "containment"
        reports intervals where the overlapping is contained within it.

    Examples
    --------
    >>> import pyranges as pr
    >>> a = '''Chromosome Start End
    ... chr1    6    12
    ... chr1    10    20
    ... chr1    22    27
    ... chr1    24    30'''

    >>> b = '''Chromosome Start End
    ... chr1    12    32
    ... chr1    14    30'''

    >>> c = '''Chromosome Start End
    ... chr1    8    15
    ... chr1    10    14
    ... chr1    32    34'''

    >>> grs = {n: pr.from_string(s) for n, s in zip(["a", "b", "c"], [a, b, c])}
    >>> for k, v in grs.items():
    ...     print("Name: " + k)
    ...     print(v)
    Name: a
    Chromosome      Start      End
    object          int64    int64
    ------------  -------  -------
    chr1                6       12
    chr1               10       20
    chr1               22       27
    chr1               24       30
    PyRanges with 4 rows and 3 columns.
    Contains 1 chromosomes.
    Name: b
    Chromosome      Start      End
    object          int64    int64
    ------------  -------  -------
    chr1               12       32
    chr1               14       30
    PyRanges with 2 rows and 3 columns.
    Contains 1 chromosomes.
    Name: c
    Chromosome      Start      End
    object          int64    int64
    ------------  -------  -------
    chr1                8       15
    chr1               10       14
    chr1               32       34
    PyRanges with 3 rows and 3 columns.
    Contains 1 chromosomes.

    >>> pr.count_overlaps(grs)
    Chromosome    Start    End      a        b        c
    object        int64    int64    int64    int64    int64
    ------------  -------  -------  -------  -------  -------
    chr1          6        8        1        0        0
    chr1          8        10       1        0        1
    chr1          10       12       2        0        2
    chr1          12       14       1        1        2
    ...           ...      ...      ...      ...      ...
    chr1          24       27       2        2        0
    chr1          27       30       1        2        0
    chr1          30       32       0        1        0
    chr1          32       34       0        0        1
    PyRanges with 12 rows and 6 columns.
    Contains 1 chromosomes.

    >>> gr = pr.PyRanges({"Chromosome": ["chr1"], "Start": [0], "End": [40]}).tile(10)
    >>> gr
    Chromosome      Start      End
    object          int64    int64
    ------------  -------  -------
    chr1                0       10
    chr1               10       20
    chr1               20       30
    chr1               30       40
    PyRanges with 4 rows and 3 columns.
    Contains 1 chromosomes.

    >>> pr.count_overlaps(grs, gr)
    Chromosome      Start      End        a        b        c
    object          int64    int64    int64    int64    int64
    ------------  -------  -------  -------  -------  -------
    chr1                0       10        5        5        4
    chr1               10       20        5        5        4
    chr1               20       30        5        5        4
    chr1               30       40        5        5        4
    PyRanges with 4 rows and 6 columns.
    Contains 1 chromosomes.
    """
    if features is None:
        features = concat.concat(list(grs.values())).split(between=True)
    else:
        features = features.copy()

    from pyranges.methods.overlap import _count_overlaps

    for name, gr in grs.items():
        _gr = gr[[c for c in gr.columns if c in GENOME_LOC_COLS_WITH_STRAND]]
        features = features.apply_pair(
            _gr, _count_overlaps, by=by, strand_behavior=strand_behavior, how=how, return_indexes=True, name=name
        )
    return features.astype({k: int for k in grs.keys()})
