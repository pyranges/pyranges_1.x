"""Module for PyRanges concat method."""

from collections.abc import Iterable
from typing import TYPE_CHECKING

import pandas as pd

import pyranges as pr
if TYPE_CHECKING:
    from pyranges import PyRanges


def concat(grs: Iterable["PyRanges"], *args, **kwargs) -> "PyRanges":
    """Concatenate PyRanges.

    Parameters
    ----------
    grs

    Returns
    -------
    pyranges.PyRanges

    Examples
    --------
    >>> gr1 = pr.example_data.f2
    >>> gr2 = pr.example_data.f1
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
    return pr.PyRanges(pd.concat(grs, *args, **kwargs))
