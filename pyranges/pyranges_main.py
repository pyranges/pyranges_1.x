"""Data structure for genomic intervals and their annotation."""
import heapq
import itertools
import shutil
import typing
from collections import defaultdict
from functools import cached_property
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Optional,
    Self,
    Tuple,
    Union,
)

import numpy as np
import pandas as pd
from natsort import natsorted  # type: ignore
from tabulate import tabulate

import pyranges as pr
from pyranges.methods.intersection import _intersection, _overlap
from pyranges.multithreaded import (
    _extend,
    _extend_grp,
    _tes,
    _tss,
    pyrange_apply,
    pyrange_apply_single,
)
from pyranges.names import (
    STRAND_COL,
    FORWARD_STRAND,
    REVERSE_STRAND,
    GENOME_LOC_COLS,
    GENOME_LOC_COLS_WITH_STRAND,
    VALID_GENOMIC_STRAND_INFO,
    VALID_STRAND_BEHAVIOR_TYPE,
    VALID_OVERLAP_TYPE, CHROM_COL,
)
from pyranges.tostring import adjust_table_width

if TYPE_CHECKING:
    from pathlib import Path

    from pandas.core.indexes.base import Index
    from pyrle.rledict import RleDict  # type: ignore

__all__ = ["PyRanges"]


ChromosomeLocation = Union[str, Tuple[str, str]]


def fill_kwargs(kwargs: Dict[str, Any]) -> Dict[str, Any]:
    """Give the kwargs dict default options."""

    defaults = {
        "strand_behavior": None,
        "overlap": True,
        "how": None,
        "invert": None,
        "new_pos": None,
        "suffixes": ["_a", "_b"],
        "suffix": "_b",
        "sparse": {"self": False, "other": False},
    }

    defaults.update(kwargs)

    return defaults


class PyRanges(pr.RangeFrame):

    """Two-dimensional representation of genomic intervals and their annotations.

    A PyRanges object must have the columns Chromosome, Start and End. These
    describe the genomic position and function as implicit row labels. A Strand
    column is optional and adds strand information to the intervals. Any other
    columns are allowed and are considered metadata.

    Operations between PyRanges align intervals based on their position.

    If a PyRanges is built using the arguments chromosomes, starts, ends and
    optionally strands, all non-scalars must be of the same length.

    Parameters
    ----------
    df : DataFrame or dict of DataFrame, default None
        The data to be stored in the PyRanges.

    chromosomes : array-like or scalar value, default None
        The chromosome(s) in the PyRanges.

    starts : array-like, default None
        The start postions in the PyRanges.

    ends : array-like, default None
        The end postions in the PyRanges.

    strands : array-like or scalar value, default None
        The strands in the PyRanges.

    copy_df : bool, default True
        Copy input DataFrame

    See Also
    --------

    pyranges.read_bed: read bed-file into PyRanges
    pyranges.read_bam: read bam-file into PyRanges
    pyranges.read_gff: read gff-file into PyRanges
    pyranges.read_gtf: read gtf-file into PyRanges
    pyranges.from_string: create PyRanges from multiline string

    Notes
    -----

    A PyRanges object is represented internally as a dictionary efficiency. The keys are
    chromosomes or chromosome/strand tuples and the values are pandas pd.DataFrames.

    Examples
    --------

    >>> pr.PyRanges()
    Empty PyRanges

    >>> pr.from_args(chromosomes="chr1", starts=(1, 5), ends=[3, 149],
    ...             strands=("+", "-"))
    +--------------+-----------+-----------+--------------+
    | Chromosome   |     Start |       End | Strand       |
    | (category)   |   (int64) |   (int64) | (category)   |
    |--------------+-----------+-----------+--------------|
    | chr1         |         1 |         3 | +            |
    | chr1         |         5 |       149 | -            |
    +--------------+-----------+-----------+--------------+
    Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    >>> df = pd.DataFrame({"Chromosome": ["chr1", "chr2"], "Start": [100, 200],
    ...                    "End": [150, 201]})
    >>> df
      Chromosome  Start  End
    0       chr1    100  150
    1       chr2    200  201
    >>> pr.PyRanges(df)
    +--------------+-----------+-----------+
    | Chromosome   |     Start |       End |
    | (category)   |   (int64) |   (int64) |
    |--------------+-----------+-----------|
    | chr1         |       100 |       150 |
    | chr2         |       200 |       201 |
    +--------------+-----------+-----------+
    Unstranded PyRanges object has 2 rows and 3 columns from 2 chromosomes.
    For printing, the PyRanges was sorted on Chromosome.


    >>> gr = pr.PyRanges({"Chromosome": [1, 1], "Strand": ["+", "-"], "Start": [1, 4], "End": [2, 27],
    ...                    "TP": [0, 1], "FP": [12, 11], "TN": [10, 9], "FN": [2, 3]})
    >>> gr
    +--------------+--------------+-----------+-----------+-----------+-----------+-----------+-----------+
    |   Chromosome | Strand       |     Start |       End |        TP |        FP |        TN |        FN |
    |   (category) | (category)   |   (int64) |   (int64) |   (int64) |   (int64) |   (int64) |   (int64) |
    |--------------+--------------+-----------+-----------+-----------+-----------+-----------+-----------|
    |            1 | +            |         1 |         2 |         0 |        12 |        10 |         2 |
    |            1 | -            |         4 |        27 |         1 |        11 |         9 |         3 |
    +--------------+--------------+-----------+-----------+-----------+-----------+-----------+-----------+
    Stranded PyRanges object has 2 rows and 8 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    dfs: Union[Dict[str, pd.DataFrame], Dict[Tuple[str, str], pd.DataFrame]]
    """Dict mapping chromosomes or chromosome/strand pairs to pandas pd.DataFrames."""

    features = None
    """Namespace for genomic-features methods.

    See Also
    --------
    pyranges.genomicfeatures : namespace for feature-functionality
    pyranges.genomicfeatures.GenomicFeaturesMethods : namespace for feature-functionality
    """

    stats = None
    """Namespace for statistcal methods.

    See Also
    --------
    pyranges.statistics : namespace for statistics
    pyranges.stats.StatisticsMethods : namespace for statistics
    """

    @cached_property
    def _required_columns(self) -> Iterable[str]:
        return GENOME_LOC_COLS[:]

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        missing_columns = set(self._required_columns) - set(self.columns)
        if missing_columns:
            msg = f"Missing required columns: {missing_columns}"
            raise ValueError(msg)
        self._check_index_column_names()

    def _chrom_and_strand_info(self) -> str:
        num_chrom = self[CHROM_COL].nunique()
        strand_info = ""

        if self.has_strand_column():
            num_strands = self[STRAND_COL].nunique()
            strands = f" and {num_strands} strands"
            if not self.strand_values_valid:
                nongenomic_strands = sorted(set(self[STRAND_COL]).difference(VALID_GENOMIC_STRAND_INFO))
                example_strands = ", ".join(
                    [str(s) for s in nongenomic_strands[:7]] + (["..."] if len(nongenomic_strands) > 7 else []))
                strands += f" (including non-genomic strands: {example_strands})"
            strand_info = strands

        return f"Contains {num_chrom} chromosomes{strand_info}"

    def __str__(self, max_col_width: int | None = None, max_total_width: int | None = None) -> str:
        """
        Return string representation.

        Examples
        --------
        >>> gr = pr.PyRanges({"Chromosome": [3, 2, 1], "Start": [0, 100, 250], "End": [10, 125, 251]})
        >>> gr
          Chromosome    Start      End
               int64    int64    int64
        ------------  -------  -------
                   3        0       10
                   2      100      125
                   1      250      251
        PyRanges with 3 rows and 3 columns.
        Contains 3 chromosomes.
        >>> gr.insert(3, "Strand", ["+", "-", "+"])
        >>> gr
          Chromosome    Start      End  Strand
               int64    int64    int64  object
        ------------  -------  -------  --------
                   3        0       10  +
                   2      100      125  -
                   1      250      251  +
        PyRanges with 3 rows and 4 columns.
        Contains 3 chromosomes and 2 strands.
        >>> gr.loc[:, "Strand"] = [".", "^", "/"]
        >>> gr
          Chromosome    Start      End  Strand
               int64    int64    int64  object
        ------------  -------  -------  --------
                   3        0       10  .
                   2      100      125  ^
                   1      250      251  /
        PyRanges with 3 rows and 4 columns.
        Contains 3 chromosomes and 3 strands (including non-genomic strands: ., /, ^).
        >>> gr = pr.concat([gr, gr, gr])
        >>> gr  # If a PyRanges has more than eight rows the repr is truncated in the middle.
        Chromosome    Start    End      Strand
        int64         int64    int64    object
        ------------  -------  -------  --------
        3             0        10       .
        2             100      125      ^
        1             250      251      /
        3             0        10       .
        ...           ...      ...      ...
        1             250      251      /
        3             0        10       .
        2             100      125      ^
        1             250      251      /
        PyRanges with 9 rows and 4 columns.
        Contains 3 chromosomes and 3 strands (including non-genomic strands: ., /, ^).
        >>> gr = PyRanges({"Chromosome": [1], "Start": [1], "End": [2], "Strand": ["+"], "Name": ["Sonic The Hedgehog"], "gene_id": ["ENSG00000188976"], "short_gene_name": ["SHH"], "type": ["transcript"]})
        >>> str_repr = gr.__str__(max_col_width=10, max_total_width=80)
        >>> print(str_repr)
          Chromosome    Start      End  Strand    Name        gene_id     ...
               int64    int64    int64  object    object      object      ...
        ------------  -------  -------  --------  ----------  ----------  -----
                   1        1        2  +         Sonic T...  ENSG000...  ...
        PyRanges with 1 rows and 8 columns (2 columns not shown: "short_gene_name", "type").
        Contains 1 chromosomes and 1 strands.
        >>> gr.insert(gr.shape[1], "Score", int(1e100))
        >>> gr.insert(gr.shape[1], "Score2", int(1e100))
        >>> str_repr = gr.__str__(max_col_width=10, max_total_width=80)
        >>> print(str_repr)
          Chromosome    Start      End  Strand    Name        gene_id     ...
               int64    int64    int64  object    object      object      ...
        ------------  -------  -------  --------  ----------  ----------  -----
                   1        1        2  +         Sonic T...  ENSG000...  ...
        PyRanges with 1 rows and 10 columns (4 columns not shown: "short_gene_name", "type", "Score", ...).
        Contains 1 chromosomes and 1 strands.
        """
        return f"{super().__str__(max_total_width=max_total_width, max_col_width=max_col_width)}\n{self._chrom_and_strand_info()}."

    def _check_index_column_names(self):
        # Check for columns with the same name as the index
        if self.index.name is not None and self.index.name in self.columns:
            raise ValueError(
                f"Index name '{self.index.name}' cannot be the same as a column name."
            )
        if isinstance(self.index, pd.MultiIndex):
            for level in self.index.names:
                if level in self.columns:
                    raise ValueError(
                        f"Index level name '{level}' cannot be the same as a column name."
                    )

    def set_index(self, *args, **kwargs) -> Optional["PyRanges"]:
        # Custom behavior for set_index
        result = PyRanges(super().set_index(*args, **kwargs))
        if kwargs.get("inplace"):
            self._check_index_column_names()
        else:
            result._check_index_column_names()
        return result

    def __array_ufunc__(self, *args, **kwargs) -> "PyRanges":
        """Apply unary numpy-function.


        Apply function to all columns which are not index, i.e. Chromosome,
        Start, End nor Strand.

        Notes
        -----

        Function must produce a vector of equal length.

        Examples
        --------

        >>> gr = pr.PyRanges({"Chromosome": [1, 2, 3], "Start": [1, 2, 3],
        ... "End": [2, 3, 4], "Score": [9, 16, 25], "Score2": [121, 144, 169],
        ... "Name": ["n1", "n2", "n3"]})
        >>> gr
        +--------------+-----------+-----------+-----------+-----------+------------+
        |   Chromosome |     Start |       End |     Score |    Score2 | Name       |
        |   (category) |   (int64) |   (int64) |   (int64) |   (int64) | (object)   |
        |--------------+-----------+-----------+-----------+-----------+------------|
        |            1 |         1 |         2 |         9 |       121 | n1         |
        |            2 |         2 |         3 |        16 |       144 | n2         |
        |            3 |         3 |         4 |        25 |       169 | n3         |
        +--------------+-----------+-----------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 6 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> np.sqrt(gr)
        +--------------+-----------+-----------+-------------+-------------+------------+
        |   Chromosome |     Start |       End |       Score |      Score2 | Name       |
        |   (category) |   (int64) |   (int64) |   (float64) |   (float64) | (object)   |
        |--------------+-----------+-----------+-------------+-------------+------------|
        |            1 |         1 |         2 |           3 |          11 | n1         |
        |            2 |         2 |         3 |           4 |          12 | n2         |
        |            3 |         3 |         4 |           5 |          13 | n3         |
        +--------------+-----------+-----------+-------------+-------------+------------+
        Unstranded PyRanges object has 3 rows and 6 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        func, call, gr = args

        columns = list(gr.columns)
        non_index = [
            c for c in columns if c not in ["Chromosome", "Start", "End", "Strand"]
        ]

        for chromosome, df in gr:
            subset = df.head(1)[non_index].select_dtypes(include=np.number).columns
            _v = getattr(func, call)(df[subset], **kwargs)
            df[subset] = _v

        return gr

    def genome_loc(self, val: Any) -> "PyRanges":
        """Fetch columns or subset on position.

        If a list is provided, the column(s) in the list is returned. This subsets on columns.

        If a numpy array is provided, it must be of type bool and the same length as the PyRanges.

        Otherwise, a subset of the rows is returned with the location info provided.

        Parameters
        ----------
        val : bool array/pd.Series, tuple, list, str or slice

            Data to fetch.

        Examples
        --------

        >>> gr = pr.data.ensembl_gtf()
        >>> list(gr.columns)
        ['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'gene_biotype', 'gene_id', 'gene_name', 'gene_source', 'gene_version', 'tag', 'transcript_biotype', 'transcript_id', 'transcript_name', 'transcript_source', 'transcript_support_level', 'transcript_version', 'exon_id', 'exon_number', 'exon_version', '(assigned', 'previous', 'protein_id', 'protein_version', 'ccds_id']

        >>> gr = gr[["Source", "Feature", "gene_id"]]
        >>> gr
        +--------------+------------+--------------+-----------+-----------+--------------+-----------------+
        | Chromosome   | Source     | Feature      | Start     | End       | Strand       | gene_id         |
        | (category)   | (object)   | (category)   | (int64)   | (int64)   | (category)   | (object)        |
        |--------------+------------+--------------+-----------+-----------+--------------+-----------------|
        | 1            | havana     | gene         | 11868     | 14409     | +            | ENSG00000223972 |
        | 1            | havana     | transcript   | 11868     | 14409     | +            | ENSG00000223972 |
        | 1            | havana     | exon         | 11868     | 12227     | +            | ENSG00000223972 |
        | 1            | havana     | exon         | 12612     | 12721     | +            | ENSG00000223972 |
        | ...          | ...        | ...          | ...       | ...       | ...          | ...             |
        | 1            | havana     | gene         | 1173055   | 1179555   | -            | ENSG00000205231 |
        | 1            | havana     | transcript   | 1173055   | 1179555   | -            | ENSG00000205231 |
        | 1            | havana     | exon         | 1179364   | 1179555   | -            | ENSG00000205231 |
        | 1            | havana     | exon         | 1173055   | 1176396   | -            | ENSG00000205231 |
        +--------------+------------+--------------+-----------+-----------+--------------+-----------------+
        Stranded PyRanges object has 2,446 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        Create boolean pd.Series and use it to subset:

        >>> s = (gr.Feature == "gene") | (gr.gene_id == "ENSG00000223972")
        >>> gr[s]
        +--------------+----------------+--------------+-----------+-----------+--------------+-----------------+
        | Chromosome   | Source         | Feature      | Start     | End       | Strand       | gene_id         |
        | (category)   | (object)       | (category)   | (int64)   | (int64)   | (category)   | (object)        |
        |--------------+----------------+--------------+-----------+-----------+--------------+-----------------|
        | 1            | havana         | gene         | 11868     | 14409     | +            | ENSG00000223972 |
        | 1            | havana         | transcript   | 11868     | 14409     | +            | ENSG00000223972 |
        | 1            | havana         | exon         | 11868     | 12227     | +            | ENSG00000223972 |
        | 1            | havana         | exon         | 12612     | 12721     | +            | ENSG00000223972 |
        | ...          | ...            | ...          | ...       | ...       | ...          | ...             |
        | 1            | havana         | gene         | 1062207   | 1063288   | -            | ENSG00000273443 |
        | 1            | ensembl_havana | gene         | 1070966   | 1074306   | -            | ENSG00000237330 |
        | 1            | ensembl_havana | gene         | 1081817   | 1116361   | -            | ENSG00000131591 |
        | 1            | havana         | gene         | 1173055   | 1179555   | -            | ENSG00000205231 |
        +--------------+----------------+--------------+-----------+-----------+--------------+-----------------+
        Stranded PyRanges object has 95 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> cs = pr.data.chipseq()
        >>> cs[10000:100000]
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr2         |     33241 |     33266 | U0         |         0 | +            |
        | chr2         |     13611 |     13636 | U0         |         0 | -            |
        | chr2         |     32620 |     32645 | U0         |         0 | -            |
        | chr3         |     87179 |     87204 | U0         |         0 | +            |
        | chr4         |     45413 |     45438 | U0         |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 5 rows and 6 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> cs["chr1", "-"]
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   | Start     | End       | Name       | Score     | Strand       |
        | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 100079649 | 100079674 | U0         | 0         | -            |
        | chr1         | 223587418 | 223587443 | U0         | 0         | -            |
        | chr1         | 202450161 | 202450186 | U0         | 0         | -            |
        | chr1         | 156338310 | 156338335 | U0         | 0         | -            |
        | ...          | ...       | ...       | ...        | ...       | ...          |
        | chr1         | 203557775 | 203557800 | U0         | 0         | -            |
        | chr1         | 28114107  | 28114132  | U0         | 0         | -            |
        | chr1         | 21622765  | 21622790  | U0         | 0         | -            |
        | chr1         | 80668132  | 80668157  | U0         | 0         | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 437 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> cs["chr5", "-", 90000:]
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   | Start     | End       | Name       | Score     | Strand       |
        | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr5         | 399682    | 399707    | U0         | 0         | -            |
        | chr5         | 1847502   | 1847527   | U0         | 0         | -            |
        | chr5         | 5247533   | 5247558   | U0         | 0         | -            |
        | chr5         | 5300394   | 5300419   | U0         | 0         | -            |
        | ...          | ...       | ...       | ...        | ...       | ...          |
        | chr5         | 178786234 | 178786259 | U0         | 0         | -            |
        | chr5         | 179268931 | 179268956 | U0         | 0         | -            |
        | chr5         | 179289594 | 179289619 | U0         | 0         | -            |
        | chr5         | 180513795 | 180513820 | U0         | 0         | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 285 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> cs["chrM"]
        Empty PyRanges
        """

        from pyranges.methods.getitem import _getitem

        return _getitem(self, val)

    def __iter__(self):
        """Iterate over the keys and values.

        See Also
        --------
        pyranges.iter : iterate over multiple PyRanges

        Examples
        --------
        >>> gr = pr.PyRanges({"Chromosome": [1, 1, 1], "Start": [0, 100, 250],
        ...                   "End": [10, 125, 251], "Strand": ["+", "+", "-"]})

        >>> for k, v in gr:
        ...     print(k)
        ...     print(v)
        ('1', '+')
          Chromosome  Start  End Strand
        0          1      0   10      +
        1          1    100  125      +
        ('1', '-')
          Chromosome  Start  End Strand
        2          1    250  251      -
        """

        return iter(self.items())

    def apply_general(
        self, f: Callable, strand: Optional[bool] = None, **kwargs
    ) -> Union[Dict[str, Any], Dict[Tuple[str, str], Any]]:
        """Apply a function to the PyRanges and return a dict of dict.

        Parameters
        ----------
        f : function
            Function to apply on each pd.DataFrame in a PyRanges

        strand : Optional[bool], default None, i.e. auto

            Whether to do operations on chromosome/strand pairs or chromosomes. If None, will use
            chromosome/strand pairs if the PyRanges is stranded.

        **kwargs
            Additional keyword arguments to pass as keyword arguments to `f`

        Returns
        -------
        PyRanges
            Result of applying f to each pd.DataFrame in the PyRanges

        See also
        --------

        pyranges.PyRanges.apply: apply a function to a PyRanges and return a PyRanges
        pyranges.PyRanges.apply_pair: apply a function to a pair of PyRanges and return a PyRanges
        pyranges.PyRanges.apply_pair_general: apply a function to a pair of PyRanges and return a dict

        Note
        ----

        This is the function used internally to carry out almost all unary PyRanges methods.

        Examples
        --------

        >>> gr = pr.PyRanges({"Chromosome": [1, 1, 2, 2], "Strand": ["+", "+", "-", "+"],
        ...                    "Start": [1, 4, 2, 9], "End": [2, 27, 13, 10]})

        >>> gr.apply_general(lambda df: len(df))
        {('1', '+'): 2, ('2', '+'): 1, ('2', '-'): 1}

        >>> gr.apply_general(lambda df: len(df), strand=False)
        {'1': 2, '2': 2}
        """

        if strand is None:
            strand = self.strand_values_valid

        kwargs.update({"strand": strand})
        kwargs.update(kwargs.get("kwargs", {}))
        kwargs = fill_kwargs(kwargs)

        return pyrange_apply_single(f, self, **kwargs)

    def apply_pair(
        self,
        other: "PyRanges",
        f: Callable,
        strand_behavior: Optional[str] = None,
        **kwargs,
    ) -> "PyRanges":
        """Apply a function to a pair of PyRanges.

        The function is applied to each chromosome or chromosome/strand pair found in at least one
        of the PyRanges.

        Parameters
        ----------
        f : function
            Row-based or associative function to apply on the pd.DataFrames.

        other : PyRanges

        strand_behavior : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        as_pyranges : bool, default False

            Whether to return as a PyRanges or dict. If `f` does not return a pd.DataFrame valid for
            PyRanges, `as_pyranges` must be False.

        **kwargs
            Additional keyword arguments to pass as keyword arguments to `f`

        Returns
        -------
        dict of lists
            Result of applying f to each partition of the pd.DataFrames in the PyRanges.

        See also
        --------

        pyranges.PyRanges.apply_pair: apply a function to a pair of PyRanges
        pyranges.PyRanges.apply_chunks: apply a row-based function to a PyRanges in parallel
        pyranges.iter: iterate over two or more PyRanges

        Note
        ----

        This is the function used internally to carry out almost all comparison functions in
        PyRanges.

        Examples
        --------

        >>> gr = pr.data.chipseq()
        >>> gr2 = pr.data.chipseq_background()

        >>> gr.apply_pair(gr2, pr.methods.intersection._intersection) # same as gr.intersect(gr2)
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 226987603 | 226987617 | U0         |         0 | +            |
        | chr8         |  38747236 |  38747251 | U0         |         0 | -            |
        | chr15        |  26105515 |  26105518 | U0         |         0 | +            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f1 = pr.data.f1()
        >>> f1
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f2 = pr.data.f2()
        >>> f2
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         1 |         2 | a          |         0 | +            |
        | chr1         |         6 |         7 | b          |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        kwargs.update({"strand_behavior": strand_behavior})
        kwargs.update(kwargs.get("kwargs", {}))
        kwargs = fill_kwargs(kwargs)

        return pyrange_apply(f, self, other, **kwargs)

    def apply_pair_general(
        self, other: "PyRanges", f: Callable, strand_behavior: None = None, **kwargs
    ) -> Union[Dict[str, Any], Dict[Tuple[str, str], Any]]:
        """Apply a function to a pair of PyRanges.

        The function is applied to each chromosome or chromosome/strand pair found in at least one
        of the PyRanges.

        Parameters
        ----------
        f : function
            Row-based or associative function to apply on the pd.DataFrames.

        other : PyRanges

        strand_behavior : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        **kwargs
            Additional keyword arguments to pass as keyword arguments to `f`

        Returns
        -------
        dict of lists
            Result of applying f to each partition of the pd.DataFrames in the PyRanges.

        See also
        --------

        pyranges.PyRanges.apply: apply a function to a pair of PyRanges
        pyranges.PyRanges.apply_general: apply a function to a PyRanges and return a dict of Any
        pyranges.PyRanges.apply_pair: apply a function to a pair of PyRanges
        pyranges.iter: iterate over two or more PyRanges

        Note
        ----

        This is the function used internally to carry out almost all comparison functions in
        PyRanges.

        Examples
        --------

        >>> f1 = pr.data.f1()
        >>> f1
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f2 = pr.data.f2()
        >>> f2
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         1 |         2 | a          |         0 | +            |
        | chr1         |         6 |         7 | b          |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f1.apply_pair_general(f2, lambda df, df2: (len(df), len(df2)))
        {('chr1', '+'): (2, 2), ('chr1', '-'): (1, 2)}
        """

        kwargs.update({"strand_behavior": strand_behavior})
        kwargs.update(kwargs.get("kwargs", {}))
        kwargs = fill_kwargs(kwargs)

        result = pyrange_apply(f, self, other, **kwargs)
        return result

    def assign(
        self, col: str, f: Callable, strand: Optional[bool] = None, **kwargs
    ) -> "PyRanges":
        """Add or replace a column.

        Does not change the original PyRanges.

        Parameters
        ----------

        col : str

            Name of column.

        f : function
            Function to create new column.

        strand : Optional[bool], default None, i.e. auto

            Whether to do operations on chromosome/strand pairs or chromosomes. If None, will use
            chromosome/strand pairs if the PyRanges is stranded.

        **kwargs
            Additional keyword arguments to pass as keyword arguments to `f`


        Returns
        -------
        PyRanges
            A copy of the PyRanges with the column inserted.

        Examples
        --------

        >>> gr = pr.PyRanges({"Chromosome": [1, 1], "Start": [1, 2], "End": [3, 5],
        ... "Name": ["a", "b"]})
        >>> gr
        +--------------+-----------+-----------+------------+
        |   Chromosome |     Start |       End | Name       |
        |   (category) |   (int64) |   (int64) | (object)   |
        |--------------+-----------+-----------+------------|
        |            1 |         1 |         3 | a          |
        |            1 |         2 |         5 | b          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.assign("Blabla", lambda df: df.Chromosome.astype(str) + "_yadayada")
        +--------------+-----------+-----------+------------+------------+
        |   Chromosome |     Start |       End | Name       | Blabla     |
        |   (category) |   (int64) |   (int64) | (object)   | (object)   |
        |--------------+-----------+-----------+------------+------------|
        |            1 |         1 |         3 | a          | 1_yadayada |
        |            1 |         2 |         5 | b          | 1_yadayada |
        +--------------+-----------+-----------+------------+------------+
        Unstranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        Note that assigning to an existing name replaces the column:

        >>> gr.assign("Name",
        ... lambda df, **kwargs: df.Start.astype(str) + kwargs["sep"] +
        ... df.Name.str.capitalize(), sep="_")
        +--------------+-----------+-----------+------------+
        |   Chromosome |     Start |       End | Name       |
        |   (category) |   (int64) |   (int64) | (object)   |
        |--------------+-----------+-----------+------------|
        |            1 |         1 |         3 | 1_A        |
        |            1 |         2 |         5 | 2_B        |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        self = self.copy()

        if strand is None:
            strand = self.strand_values_valid

        kwargs["strand"] = strand
        kwargs = fill_kwargs(kwargs)

        result = pyrange_apply_single(f, self, **kwargs)

        first_result = next(iter(result.values()))

        assert isinstance(
            first_result, pd.Series
        ), "result of assign function must be pd.Series, but is {}".format(
            type(first_result)
        )

        # do a deepcopy of object
        new_self = self.copy()
        new_self.__setattr__(col, result)

        return new_self

    def boundaries(
        self, group_by: str, agg: Optional[Dict[str, Union[str, Callable]]] = None
    ) -> "PyRanges":
        """Return the boundaries of groups of intervals (e.g. transcripts)

        Parameters
        ----------

        group_by : str or list of str

            Name(s) of column(s) to group intervals

        agg : dict or None

            Defines how to aggregate metadata columns. Provided as
            dictionary of column names -> functions, function names or list of such,
            as accepted by the pd.DataFrame.agg method.


        Returns
        -------
        PyRanges
            One interval per group, with the min(Start) and max(End) of the group


        Examples
        --------

        >>> d = {"Chromosome": [1, 1, 1], "Start": [1, 60, 110], "End": [40, 68, 130], "transcript_id": ["tr1", "tr1", "tr2"], "meta": ["a", "b", "c"]}
        >>> gr = pr.PyRanges(d)
        >>> gr.length=gr.lengths()
        >>> gr
        +--------------+-----------+-----------+-----------------+------------+-----------+
        |   Chromosome |     Start |       End | transcript_id   | meta       |    length |
        |   (category) |   (int64) |   (int64) | (object)        | (object)   |   (int64) |
        |--------------+-----------+-----------+-----------------+------------+-----------|
        |            1 |         1 |        40 | tr1             | a          |        39 |
        |            1 |        60 |        68 | tr1             | b          |         8 |
        |            1 |       110 |       130 | tr2             | c          |        20 |
        +--------------+-----------+-----------+-----------------+------------+-----------+
        Unstranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.boundaries("transcript_id")
        +--------------+-----------+-----------+-----------------+
        |   Chromosome |     Start |       End | transcript_id   |
        |   (category) |   (int64) |   (int64) | (object)        |
        |--------------+-----------+-----------+-----------------|
        |            1 |         1 |        68 | tr1             |
        |            1 |       110 |       130 | tr2             |
        +--------------+-----------+-----------+-----------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.boundaries("transcript_id", agg={"length":"sum", "meta": ",".join})
        +--------------+-----------+-----------+-----------------+------------+-----------+
        |   Chromosome |     Start |       End | transcript_id   | meta       |    length |
        |   (category) |   (int64) |   (int64) | (object)        | (object)   |   (int64) |
        |--------------+-----------+-----------+-----------------+------------+-----------|
        |            1 |         1 |        68 | tr1             | a,b        |        47 |
        |            1 |       110 |       130 | tr2             | c          |        20 |
        +--------------+-----------+-----------+-----------------+------------+-----------+
        Unstranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        """
        from pyranges.methods.boundaries import _bounds

        kwargs = {"group_by": group_by, "agg": agg, "strand": self.strand_values_valid}
        kwargs = fill_kwargs(kwargs)

        result = pyrange_apply_single(_bounds, self, **kwargs)
        return pr.PyRanges(result)

    def calculate_frame(self, by: Union[str, List[str]]) -> "PyRanges":
        """Calculate the frame of each genomic interval, assuming all are coding sequences (CDS), and add it as column inplace.

        After this, the input Pyranges will contain an added "Frame" column, which determines the base of the CDS that is the first base of a codon.
        Resulting values are in range between 0 and 2 included. 0 indicates that the first base of the CDS is the first base of a codon,
        1 indicates the second base and 2 indicates the third base of the CDS.
        While the 5'-most interval of each transcript has always 0 frame, the following ones may have any of these values.

        Parameters
        ----------
        by : str or list of str

            Column(s) to group by the intervals: coding exons belonging to the same transcript have the same values in this/these column(s).

        Returns
        -------
        PyRanges

        Examples
        --------
        >>> p = pr.PyRanges({"Chromosome": [1,1,1,2,2],
        ...                   "Strand": ["+","+","+","-","-"],
        ...                   "Start": [1,31,52,101,201],
        ...                   "End": [10,45,90,130,218],
        ...                   "transcript_id": ["t1","t1","t1","t2","t2"]})
        >>> p
        +--------------+--------------+-----------+-----------+-----------------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |
        |--------------+--------------+-----------+-----------+-----------------|
        |            1 | +            |         1 |        10 | t1              |
        |            1 | +            |        31 |        45 | t1              |
        |            1 | +            |        52 |        90 | t1              |
        |            2 | -            |       101 |       130 | t2              |
        |            2 | -            |       201 |       218 | t2              |
        +--------------+--------------+-----------+-----------+-----------------+
        Stranded PyRanges object has 5 rows and 5 columns from 2 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> p.calculate_frame(by=['transcript_id'])
        +--------------+--------------+-----------+-----------+-----------------+-----------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |     Frame |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |   (int64) |
        |--------------+--------------+-----------+-----------+-----------------+-----------|
        |            1 | +            |         1 |        10 | t1              |         0 |
        |            1 | +            |        31 |        45 | t1              |         9 |
        |            1 | +            |        52 |        90 | t1              |        23 |
        |            2 | -            |       101 |       130 | t2              |        17 |
        |            2 | -            |       201 |       218 | t2              |         0 |
        +--------------+--------------+-----------+-----------+-----------------+-----------+
        Stranded PyRanges object has 5 rows and 6 columns from 2 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        """
        _self = self.copy()
        # Column to save the initial index
        _self.__index__ = np.arange(len(self))

        # Filtering for desired columns
        if isinstance(by, str):
            lst = [by]
        else:
            lst = by
        sorted_p = _self[["Strand", "__index__"] + lst]

        # Sorting by 5' (Intervals on + are sorted by ascending order and - are sorted by descending order)
        sorted_p = sorted_p.sort(by="5")

        # Creating a column saving the length for the intervals (for selenoprofiles and ensembl)
        sorted_p.__length__ = sorted_p.lengths()

        # Creating a column saving the cumulative length for the intervals
        sorted_p.loc["__cumsum__", :] = sorted_p.groupby(by=by).__length__.cumsum()

        # Creating a frame column
        sorted_p.Frame = sorted_p.__cumsum__ - sorted_p.__length__

        # Appending the Frame of sorted_p by the index of p
        sorted_p = sorted_p.apply(lambda df: df.sort_values(by="__index__"))

        _self.Frame = sorted_p.Frame

        # Drop __index__ column
        return _self.apply(lambda df: df.drop("__index__", axis=1))

    @property
    def chromosomes(self) -> List[str]:
        """Return chromosomes in natsorted order."""

        if self.strand_values_valid:
            return natsorted(set([k[0] for k in self.keys()]))
        else:
            return natsorted(set([k for k in self.keys()]))

    @property
    def chromosomes_and_strands(self) -> List[Tuple[str, str]]:
        """Return chromosomes and strands in natsorted order."""

        if not self.strand_values_valid:
            raise ValueError("PyRanges is not stranded.")
        else:
            return natsorted(set(self.keys()))

    def cluster(
        self,
        strand: Optional[bool] = None,
        by: Optional[Union[List[str], str]] = None,
        slack: int = 0,
        count: bool = False,
    ) -> "PyRanges":
        """Give overlapping intervals a common id.

        Parameters
        ----------
        strand : bool, default None, i.e. auto

            Whether to ignore strand information if PyRanges is stranded.

        by : str or list, default None

            Only intervals with an equal value in column(s) `by` are clustered.

        slack : int, default 0

            Consider intervals separated by less than `slack` to be in the same cluster. If `slack`
            is negative, intervals overlapping less than `slack` are not considered to be in the
            same cluster.

        Returns
        -------
        PyRanges
            PyRanges with an ID-column "Cluster" added.

        Warning
        -------

        Bookended intervals (i.e. the End of a PyRanges interval is the Start of
        another one) are by default considered to overlap.
        Avoid this with slack=-1.

        See also
        --------

        PyRanges.merge: combine overlapping intervals into one

        Examples
        --------

        >>> gr = pr.PyRanges({"Chromosome": [1, 1, 1, 1], "Start": [1, 2, 3, 9],
        ...                    "End": [3, 3, 10, 12], "Gene": [1, 2, 3, 3]})
        >>> gr
        +--------------+-----------+-----------+-----------+
        |   Chromosome |     Start |       End |      Gene |
        |   (category) |   (int64) |   (int64) |   (int64) |
        |--------------+-----------+-----------+-----------|
        |            1 |         1 |         3 |         1 |
        |            1 |         2 |         3 |         2 |
        |            1 |         3 |        10 |         3 |
        |            1 |         9 |        12 |         3 |
        +--------------+-----------+-----------+-----------+
        Unstranded PyRanges object has 4 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.cluster()
        +--------------+-----------+-----------+-----------+-----------+
        |   Chromosome |     Start |       End |      Gene |   Cluster |
        |   (category) |   (int64) |   (int64) |   (int64) |   (int64) |
        |--------------+-----------+-----------+-----------+-----------|
        |            1 |         1 |         3 |         1 |         1 |
        |            1 |         2 |         3 |         2 |         1 |
        |            1 |         3 |        10 |         3 |         1 |
        |            1 |         9 |        12 |         3 |         1 |
        +--------------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 4 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.cluster(by="Gene", count=True)
        +--------------+-----------+-----------+-----------+-----------+-----------+
        |   Chromosome |     Start |       End |      Gene |   Cluster |     Count |
        |   (category) |   (int64) |   (int64) |   (int64) |   (int64) |   (int64) |
        |--------------+-----------+-----------+-----------+-----------+-----------|
        |            1 |         1 |         3 |         1 |         1 |         1 |
        |            1 |         2 |         3 |         2 |         2 |         1 |
        |            1 |         3 |        10 |         3 |         3 |         2 |
        |            1 |         9 |        12 |         3 |         3 |         2 |
        +--------------+-----------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 4 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        Avoid clustering bookended intervals with slack=-1:

        >>> gr.cluster(slack=-1)
        +--------------+-----------+-----------+-----------+-----------+
        |   Chromosome |     Start |       End |      Gene |   Cluster |
        |   (category) |   (int64) |   (int64) |   (int64) |   (int64) |
        |--------------+-----------+-----------+-----------+-----------|
        |            1 |         1 |         3 |         1 |         1 |
        |            1 |         2 |         3 |         2 |         1 |
        |            1 |         3 |        10 |         3 |         2 |
        |            1 |         9 |        12 |         3 |         2 |
        +--------------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 4 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr2 = pr.data.ensembl_gtf()[["Feature", "Source"]]
        >>> gr2.cluster(by=["Feature", "Source"])
        +--------------+--------------+---------------+-----------+-----------+--------------+-----------+
        | Chromosome   | Feature      | Source        | Start     | End       | Strand       | Cluster   |
        | (category)   | (category)   | (object)      | (int64)   | (int64)   | (category)   | (int64)   |
        |--------------+--------------+---------------+-----------+-----------+--------------+-----------|
        | 1            | CDS          | ensembl       | 69090     | 70005     | +            | 1         |
        | 1            | CDS          | ensembl       | 925941    | 926013    | +            | 2         |
        | 1            | CDS          | ensembl       | 925941    | 926013    | +            | 2         |
        | 1            | CDS          | ensembl       | 925941    | 926013    | +            | 2         |
        | ...          | ...          | ...           | ...       | ...       | ...          | ...       |
        | 1            | transcript   | havana_tagene | 167128    | 169240    | -            | 1142      |
        | 1            | transcript   | mirbase       | 17368     | 17436     | -            | 1143      |
        | 1            | transcript   | mirbase       | 187890    | 187958    | -            | 1144      |
        | 1            | transcript   | mirbase       | 632324    | 632413    | -            | 1145      |
        +--------------+--------------+---------------+-----------+-----------+--------------+-----------+
        Stranded PyRanges object has 2,446 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """
        _self = self.copy()
        if strand is None:
            strand = _self.strand_values_valid

        kwargs = {"strand": strand, "slack": slack, "count": count, "by": by}
        kwargs = fill_kwargs(kwargs)

        _stranded = _self.strand_values_valid
        if not strand and _stranded:
            _self.__Strand__ = _self.Strand
            _self = _self.remove_strand()

        if not by:
            from pyranges.methods.cluster import _cluster

            df = pyrange_apply_single(_cluster, _self, **kwargs)
        else:
            from pyranges.methods.cluster import _cluster_by

            kwargs["by"] = by
            df = pyrange_apply_single(_cluster_by, _self, **kwargs)

        gr = pr.PyRanges(df)

        # each chromosome got overlapping ids (0 to len). Need to make unique!
        new_dfs = {}
        first = True
        max_id = 0
        for k, v in gr.items():
            if first:
                max_id = v.Cluster.max()
                new_dfs[k] = v
                first = False
                continue

            v.loc[:, "Cluster"] += max_id
            max_id = v.Cluster.max()
            new_dfs[k] = v

        if not strand and _stranded:
            renamed = [
                d.rename(columns={"__Strand__": "Strand"}) for d in new_dfs.values()
            ]
            return PyRanges._zip_locationkey_and_data(
                new_dfs.keys(), renamed, strand=True
            )
        else:
            return PyRanges._zip_locationkey_and_data(
                new_dfs.keys(), new_dfs.values(), strand=strand
            )

    def count_overlaps(
        self,
        other: "PyRanges",
        strand_behavior: None = None,
        keep_nonoverlapping: bool = True,
        overlap_col: str = "NumberOverlaps",
    ) -> "PyRanges":
        """Count number of overlaps per interval.

        Count how many intervals in self overlap with those in other.

        Parameters
        ----------
        strand_behavior : {"same", "opposite", None, False}, default None, i.e. auto

            Whether to perform the operation on the same, opposite or no strand. Use False to
            ignore the strand. None means use "same" if both PyRanges are stranded, otherwise
            ignore.

        keep_nonoverlapping : bool, default True

            Keep intervals without overlaps.

        overlap_col : str, default "NumberOverlaps"

            Name of column with overlap counts.

        Returns
        -------
        PyRanges
            PyRanges with a column of overlaps added.

        See also
        --------

        PyRanges.coverage: find coverage of PyRanges
        pyranges.count_overlaps: count overlaps from multiple PyRanges

        Examples
        --------
        >>> f1 = pr.data.f1().drop()
        >>> f1
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         6 | +            |
        | chr1         |         8 |         9 | +            |
        | chr1         |         5 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        >>> f2 = pr.data.f2().drop()
        >>> f2
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         2 | +            |
        | chr1         |         6 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f1.count_overlaps(f2, overlap_col="Count")
        +--------------+-----------+-----------+--------------+-----------+
        | Chromosome   |     Start |       End | Strand       |     Count |
        | (category)   |   (int64) |   (int64) | (category)   |   (int64) |
        |--------------+-----------+-----------+--------------+-----------|
        | chr1         |         3 |         6 | +            |         0 |
        | chr1         |         8 |         9 | +            |         0 |
        | chr1         |         5 |         7 | -            |         1 |
        +--------------+-----------+-----------+--------------+-----------+
        Stranded PyRanges object has 3 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        kwargs = {
            "strand_behavior": strand_behavior,
            "keep_nonoverlapping": keep_nonoverlapping,
            "overlap_col": overlap_col,
        }
        kwargs = fill_kwargs(kwargs)

        from pyranges.methods.coverage import _number_overlapping

        counts = pyrange_apply(_number_overlapping, self, other, **kwargs)

        return pr.PyRanges(counts)

    def coverage(
        self,
        other: "PyRanges",
        strand_behavior: None = None,
        keep_nonoverlapping: bool = True,
        overlap_col: str = "NumberOverlaps",
        fraction_col: str = "FractionOverlaps",
    ) -> "PyRanges":
        """Count number of overlaps and their fraction per interval.

        Count how many intervals in self overlap with those in other.

        Parameters
        ----------
        strand_behavior : {"same", "opposite", None, False}, default None, i.e. auto

            Whether to perform the operation on the same, opposite or no strand. Use False to
            ignore the strand. None means use "same" if both PyRanges are stranded, otherwise
            ignore.

        keep_nonoverlapping : bool, default True

            Keep intervals without overlaps.

        overlap_col : str, default "NumberOverlaps"

            Name of column with overlap counts.

        fraction_col : str, default "FractionOverlaps"

            Name of column with fraction of counts.

        Returns
        -------
        PyRanges
            PyRanges with a column of overlaps added.

        See also
        --------

        pyranges.count_overlaps: count overlaps from multiple PyRanges

        Examples
        --------
        >>> f1 = pr.PyRanges({"Chromosome": [1, 1, 1], "Start": [3, 8, 5],
        ...                    "End": [6,  9, 7]})
        >>> f1
        +--------------+-----------+-----------+
        |   Chromosome |     Start |       End |
        |   (category) |   (int64) |   (int64) |
        |--------------+-----------+-----------|
        |            1 |         3 |         6 |
        |            1 |         8 |         9 |
        |            1 |         5 |         7 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        >>> f2 = pr.PyRanges({"Chromosome": [1, 1], "Start": [1, 6],
        ...                    "End": [2, 7]})
        >>> f2
        +--------------+-----------+-----------+
        |   Chromosome |     Start |       End |
        |   (category) |   (int64) |   (int64) |
        |--------------+-----------+-----------|
        |            1 |         1 |         2 |
        |            1 |         6 |         7 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> f1.coverage(f2, overlap_col="C", fraction_col="F")
        +--------------+-----------+-----------+-----------+-------------+
        |   Chromosome |     Start |       End |         C |           F |
        |   (category) |   (int64) |   (int64) |   (int64) |   (float64) |
        |--------------+-----------+-----------+-----------+-------------|
        |            1 |         3 |         6 |         0 |         0   |
        |            1 |         8 |         9 |         0 |         0   |
        |            1 |         5 |         7 |         1 |         0.5 |
        +--------------+-----------+-----------+-----------+-------------+
        Unstranded PyRanges object has 3 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        kwargs = {
            "strand_behavior": strand_behavior,
            "keep_nonoverlapping": keep_nonoverlapping,
            "overlap_col": overlap_col,
            "fraction_col": fraction_col,
        }
        kwargs = fill_kwargs(kwargs)

        counts = self.count_overlaps(
            other,
            keep_nonoverlapping=True,
            overlap_col=overlap_col,
            strand_behavior=strand_behavior,
        )

        strand = True if kwargs["strand_behavior"] else False
        other = other.merge_overlaps(count=True, strand=strand)

        from pyranges.methods.coverage import _coverage

        counts = pr.PyRanges(pyrange_apply(_coverage, counts, other, **kwargs))

        return counts

    def drop(self, *args, **kwargs) -> "PyRanges | None":
        res = super().drop(*args, **kwargs)
        if res is not None:
            return PyRanges(res)

    def drop_duplicate_positions(
        self, strand: Optional[bool] = None, keep: Union[bool, str] = "first"
    ) -> "PyRanges":
        """Return PyRanges with duplicate postion rows removed.

        Parameters
        ----------

        strand : bool, default None, i.e. auto

            Whether to take strand-information into account when considering duplicates.

        keep : {"first", "last", False}

            Whether to keep first, last or drop all duplicates.

        Examples
        --------

        >>> gr = pr.from_string('''Chromosome Start End Strand Name
        ... 1 1 2 + A
        ... 1 1 2 - B
        ... 1 1 2 + Z''')
        >>> gr
        +--------------+-----------+-----------+--------------+------------+
        |   Chromosome |     Start |       End | Strand       | Name       |
        |   (category) |   (int64) |   (int64) | (category)   | (object)   |
        |--------------+-----------+-----------+--------------+------------|
        |            1 |         1 |         2 | +            | A          |
        |            1 |         1 |         2 | +            | Z          |
        |            1 |         1 |         2 | -            | B          |
        +--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 3 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.drop_duplicate_positions()
        +--------------+-----------+-----------+--------------+------------+
        |   Chromosome |     Start |       End | Strand       | Name       |
        |   (category) |   (int64) |   (int64) | (category)   | (object)   |
        |--------------+-----------+-----------+--------------+------------|
        |            1 |         1 |         2 | +            | A          |
        |            1 |         1 |         2 | -            | B          |
        +--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.drop_duplicate_positions(keep="last")
        +--------------+-----------+-----------+--------------+------------+
        |   Chromosome |     Start |       End | Strand       | Name       |
        |   (category) |   (int64) |   (int64) | (category)   | (object)   |
        |--------------+-----------+-----------+--------------+------------|
        |            1 |         1 |         2 | +            | Z          |
        |            1 |         1 |         2 | -            | B          |
        +--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        Note that the reverse strand is considered to be behind the forward strand:

        >>> gr.drop_duplicate_positions(keep="last", strand=False)
        +--------------+-----------+-----------+--------------+------------+
        |   Chromosome |     Start |       End | Strand       | Name       |
        |   (category) |   (int64) |   (int64) | (category)   | (object)   |
        |--------------+-----------+-----------+--------------+------------|
        |            1 |         1 |         2 | -            | B          |
        +--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 1 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.drop_duplicate_positions(keep=False, strand=False)
        Empty PyRanges
        """

        from pyranges.methods.drop_duplicates import _drop_duplicate_positions

        if strand is None:
            strand = self.strand_values_valid

        kwargs = {
            "sparse": {"self": False},
            "keep": keep,
            "strand": strand and self.strand_values_valid,
        }
        kwargs = fill_kwargs(kwargs)
        return pr.PyRanges(
            pyrange_apply_single(_drop_duplicate_positions, self, **kwargs)
        )

    def extend(
        self, ext: Union[Dict[str, int], int], group_by: None = None
    ) -> "PyRanges":
        """Extend the intervals from the ends.

        Parameters
        ----------

        ext : int or dict of ints with "3" and/or "5" as keys.

            The number of nucleotides to extend the ends with.
            If an int is provided, the same extension is applied to both
            the start and end of intervals, while a dict input allows to control
            differently the two ends. Note also that 5' and 3' extensions take
            the strand into account, if the intervals are stranded.

        group_by : str or list of str, default: None

            group intervals by these column name(s), so that the extension is applied
            only to the left-most and/or right-most interval.

        See Also
        --------
        PyRanges.subsequence : obtain subsequences of intervals
        PyRanges.spliced_subsequence : obtain subsequences of intervals, providing transcript-level coordinates

        Examples
        --------

        >>> d = {'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 8, 5], 'End': [6, 9, 7],
        ...      'Strand': ['+', '+', '-']}
        >>> gr = pr.PyRanges(d)
        >>> gr
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         6 | +            |
        | chr1         |         8 |         9 | +            |
        | chr1         |         5 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.extend(4)
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         0 |        10 | +            |
        | chr1         |         4 |        13 | +            |
        | chr1         |         1 |        11 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.extend({"3": 1})
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         7 | +            |
        | chr1         |         8 |        10 | +            |
        | chr1         |         4 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.extend({"3": 1, "5": 2})
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         7 | +            |
        | chr1         |         6 |        10 | +            |
        | chr1         |         4 |         9 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.extend(-1)
        Traceback (most recent call last):
        ...
        AssertionError: Some intervals are negative or zero length after applying extend!
        """

        if isinstance(ext, dict):
            assert (
                self.strand_values_valid
            ), "PyRanges must be stranded to add 5/3-end specific extend."

        kwargs = fill_kwargs(
            {"ext": ext, "strand": self.strand_values_valid, "group_by": group_by}
        )
        func = _extend if group_by is None else _extend_grp
        return PyRanges(pyrange_apply_single(func, self, **kwargs))

    # # TODO: use subtract code here instead, easier
    # def no_overlap(self, other, **kwargs):

    #     kwargs = fill_kwargs(kwargs)
    #     kwargs["invert"] = True
    #     kwargs["sparse"] = {"self": False, "other": True}

    #     # if kwargs["strand_behavior"] in ["same", "opposite"]:
    #     #     kwargs["strand_behavior"] = {
    #     #         "same": "opposite",
    #     #         "opposite": "same"
    #     #     }[kwargs["strand_behavior"]]
    #     dfs = pyrange_apply(_overlap, self, other, **kwargs)

    #     return PyRanges(dfs)

    # @profile

    def five_end(self) -> "PyRanges":
        """Return the five prime end of intervals.

        The five prime end is the start of a forward strand or the end of a reverse strand.

        Returns
        -------
        PyRanges

            PyRanges with the five prime ends

        Notes
        -----

        Requires the PyRanges to be stranded.

        See Also
        --------

        PyRanges.three_end : return the 3' end

        Examples
        --------

        >>> gr = pr.PyRanges({'Chromosome': ['chr1', 'chr1'], 'Start': [3, 5], 'End': [9, 7],
        ...                    'Strand': ["+", "-"]})
        >>> gr
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         9 | +            |
        | chr1         |         5 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.five_end()
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         4 | +            |
        | chr1         |         6 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        assert self.strand_values_valid, "Need stranded pyrange to find 5'."
        kwargs = fill_kwargs({"strand": self.strand_values_valid})
        return pr.PyRanges(pyrange_apply_single(_tss, self, **kwargs))

    def has_strand_column(self) -> bool:
        return STRAND_COL in self.columns

    def head(self, n: int = 8) -> "PyRanges":
        """Return the n first rows.

        Parameters
        ----------

        n : int, default 8

            Return n rows.

        Returns
        -------
        PyRanges

            PyRanges with the n first rows.

        See Also
        --------

        PyRanges.tail : return the last rows
        PyRanges.sample : return random rows

        Examples
        --------

        >>> gr = pr.data.chipseq()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   | Start     | End       | Name       | Score     | Strand       |
        | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 212609534 | 212609559 | U0         | 0         | +            |
        | chr1         | 169887529 | 169887554 | U0         | 0         | +            |
        | chr1         | 216711011 | 216711036 | U0         | 0         | +            |
        | chr1         | 144227079 | 144227104 | U0         | 0         | +            |
        | ...          | ...       | ...       | ...        | ...       | ...          |
        | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
        | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
        | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
        | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.head(3)
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 212609534 | 212609559 | U0         |         0 | +            |
        | chr1         | 169887529 | 169887554 | U0         |         0 | +            |
        | chr1         | 216711011 | 216711036 | U0         |         0 | +            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        subsetter = np.zeros(len(self), dtype=np.bool_)
        subsetter[:n] = True
        return self[subsetter]

    def intersect(
        self,
        other: "PyRanges",
        strand_behavior: Optional[bool] = None,
        how: Optional[str] = None,
        invert: bool = False,
    ) -> "PyRanges":
        """Return overlapping subintervals.

        Returns the segments of the intervals in self which overlap with those in other.

        Parameters
        ----------
        other : PyRanges

            PyRanges to intersect.

        strand_behavior : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        how : {None, "first", "last", "containment"}, default None, i.e. all

            What intervals to report. By default reports all overlapping intervals. "containment"
            reports intervals where the overlapping is contained within it.

        invert : bool, default False

            Whether to return the intervals without overlaps.

        Returns
        -------
        PyRanges

            A PyRanges with overlapping subintervals.

        See also
        --------

        PyRanges.set_intersect : set-intersect PyRanges
        PyRanges.overlap : report overlapping intervals

        Examples
        --------

        >>> gr = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [1, 4, 10],
        ...                    "End": [3, 9, 11], "ID": ["a", "b", "c"]})
        >>> gr
        Chromosome      Start      End  ID
        object          int64    int64  object
        ------------  -------  -------  --------
        chr1                1        3  a
        chr1                4        9  b
        chr1               10       11  c


        >>> gr2 = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [2, 2, 9], "End": [3, 9, 10]})
        >>> gr2
        Chromosome      Start      End
        object          int64    int64
        ------------  -------  -------
        chr1                2        3
        chr1                2        9
        chr1                9       10

        >>> gr.intersect(gr2)
        Chromosome      Start      End  ID
        object          int64    int64  object
        ------------  -------  -------  --------
        chr1                2        3  a
        chr1                2        3  a
        chr1                4        9  b

        >>> gr.intersect(gr2, how="first")
        Chromosome      Start      End  ID
        object          int64    int64  object
        ------------  -------  -------  --------
        chr1                2        3  a
        chr1                4        9  b

        >>> gr.intersect(gr2, how="containment")
        Chromosome      Start      End  ID
        object          int64    int64  object
        ------------  -------  -------  --------
        chr1                4        9  b
        """

        kwargs = {
            "how": how,
            "strand_behavior": strand_behavior,
            "sparse": {"self": False, "other": True},
        }
        kwargs = fill_kwargs(kwargs)

        if len(self) == 0:
            return self

        if invert:
            self.__ix__ = np.arange(len(self))

        df = pyrange_apply(_intersection, self, other, **kwargs)
        result = pr.PyRanges(df)

        if invert:
            found_idxs = getattr(result, "__ix__", [])
            _result = self[~pd.Series(self.__ix__).isin(found_idxs)]
            return PyRanges(_result.drop("__ix__"))

        return result

    def join_overlaps(
        self,
        other: "PyRanges",
        strand_behavior: Optional[str] = None,
        how: Optional[str] = None,
        report_overlap: bool = False,
        slack: int = 0,
        suffix: str = "_b",
        apply_strand_suffix: None = None,
        preserve_order: bool = False,
    ) -> "PyRanges":
        """Join PyRanges on genomic location.

        Parameters
        ----------
        other : PyRanges

            PyRanges to join.

        strand_behavior : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        how : {None, "left", "right"}, default None, i.e. "inner"

            How to handle intervals without overlap. None means only keep overlapping intervals.
            "left" keeps all intervals in self, "right" keeps all intervals in other.

        report_overlap : bool, default False

            Report amount of overlap in base pairs.

        slack : int, default 0

            Lengthen intervals in self before joining.

        suffix : str or tuple, default "_b"

            Suffix to give overlapping columns in other.

        apply_strand_suffix : bool, default None

            If first pyranges is unstranded, but the second is not, the first will be given a strand column.
            apply_strand_suffix makes the added strand column a regular data column instead by adding a suffix.

        preserve_order : bool, default False

            If True, preserves the order after performing the join (only relevant in "outer", "left" and "right" joins).

        Returns
        -------
        PyRanges

            A PyRanges appended with columns of another.

        Notes
        -----

        The chromosome from other will never be reported as it is always the same as in self.

        As pandas did not have NaN for non-float datatypes until recently, "left" and "right" join
        give non-overlapping rows the value -1 to avoid promoting columns to object. This will
        change to NaN in a future version as general NaN becomes stable in pandas.

        See also
        --------

        PyRanges.new_position : give joined PyRanges new coordinates

        Examples
        --------

        >>> f1 = pr.PyRanges({'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 8, 5],
        ...                    'End': [6, 9, 7], 'Name': ['interval1', 'interval3', 'interval2']})
        >>> f1
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | Name       |
        | (category)   |   (int64) |   (int64) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         3 |         6 | interval1  |
        | chr1         |         8 |         9 | interval3  |
        | chr1         |         5 |         7 | interval2  |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> f2 = pr.PyRanges({'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...                    'End': [2, 7], 'Name': ['a', 'b']})
        >>> f2
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | Name       |
        | (category)   |   (int64) |   (int64) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         1 |         2 | a          |
        | chr1         |         6 |         7 | b          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> f1.join_overlaps(f2)
        +--------------+-----------+-----------+------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | Name       |   Start_b |     End_b | Name_b     |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) |   (int64) | (object)   |
        |--------------+-----------+-----------+------------+-----------+-----------+------------|
        | chr1         |         5 |         7 | interval2  |         6 |         7 | b          |
        +--------------+-----------+-----------+------------+-----------+-----------+------------+
        Unstranded PyRanges object has 1 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> f1.join_overlaps(f2, how="right")
        +--------------+-----------+-----------+------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | Name       |   Start_b |     End_b | Name_b     |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) |   (int64) | (object)   |
        |--------------+-----------+-----------+------------+-----------+-----------+------------|
        | chr1         |         5 |         7 | interval2  |         6 |         7 | b          |
        | chr1         |        -1 |        -1 | -1         |         1 |         2 | a          |
        +--------------+-----------+-----------+------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        With slack 1, bookended features are joined (see row 1):

        >>> f1.join_overlaps(f2, slack=1)
        +--------------+-----------+-----------+------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | Name       |   Start_b |     End_b | Name_b     |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) |   (int64) | (object)   |
        |--------------+-----------+-----------+------------+-----------+-----------+------------|
        | chr1         |         3 |         6 | interval1  |         6 |         7 | b          |
        | chr1         |         5 |         7 | interval2  |         6 |         7 | b          |
        +--------------+-----------+-----------+------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> f1.join_overlaps(f2, how="right", preserve_order=True)
        +--------------+-----------+-----------+------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | Name       |   Start_b |     End_b | Name_b     |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) |   (int64) | (object)   |
        |--------------+-----------+-----------+------------+-----------+-----------+------------|
        | chr1         |        -1 |        -1 | -1         |         1 |         2 | a          |
        | chr1         |         5 |         7 | interval2  |         6 |         7 | b          |
        +--------------+-----------+-----------+------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        from pyranges.methods.join import _both_dfs

        kwargs: Dict[str, Any] = {
            "strand_behavior": strand_behavior,
            "how": how,
            "report_overlap": report_overlap,
            "suffix": suffix,
            "apply_strand_suffix": apply_strand_suffix,
            "preserve_order": preserve_order,
        }
        if slack:
            self = self.copy()
            self.insert(len(self.columns), "Start__slack", self.Start)
            self.insert(len(self.columns), "End__slack", self.End)

            self = self.extend(slack)

        if "suffix" in kwargs and isinstance(kwargs["suffix"], str):
            suffixes = "", kwargs["suffix"]
            kwargs["suffixes"] = suffixes

        kwargs = fill_kwargs(kwargs)

        gr = PyRanges(pyrange_apply(_both_dfs, self, other, **kwargs))

        if slack and len(gr) > 0:
            gr.Start = gr.Start__slack
            gr.End = gr.End__slack
            gr = pr.PyRanges(gr.drop(["Start__slack", "End__slack"], axis=1))

        return gr

    @property
    def length(self) -> int:
        """Return the total length of the intervals.

        See Also
        --------

        PyRanges.lengths : return the intervals lengths

        Examples
        --------

        >>> gr = pr.data.f1()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.length
        6

        To find the length of the genome covered by the intervals, use merge first:

        >>> gr.merge_overlaps(strand=False).length
        5
        """

        lengths = self.lengths()
        length = lengths.sum()
        return int(length)

    def lengths(self) -> pd.Series:
        """Return the length of each interval.

        Parameters
        ----------

        Returns
        -------
        pd.Series or dict of pd.Series with the lengths of each interval.

        See Also
        --------

        PyRanges.lengths : return the intervals lengths

        Examples
        --------

        >>> gr = pr.data.f1()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.lengths()
        0    3
        1    1
        2    2
        dtype: int64

        To find the length of the genome covered by the intervals, use merge first:

        >>> gr.Length = gr.lengths()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+-----------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |    Length |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |   (int64) |
        |--------------+-----------+-----------+------------+-----------+--------------+-----------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |         3 |
        | chr1         |         8 |         9 | interval3  |         0 | +            |         1 |
        | chr1         |         5 |         7 | interval2  |         0 | -            |         2 |
        +--------------+-----------+-----------+------------+-----------+--------------+-----------+
        Stranded PyRanges object has 3 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        return self.End - self.Start

    def max_disjoint(
        self, strand: Optional[bool] = None, slack: int = 0, **kwargs
    ) -> "PyRanges":
        """Find the maximal disjoint set of intervals.

        Parameters
        ----------
        strand : bool, default None, i.e. auto

            Find the max disjoint set separately for each strand.

        slack : int, default 0

            Consider intervals within a distance of slack to be overlapping.

        Returns
        -------
        PyRanges

            PyRanges with maximal disjoint set of intervals.

        Examples
        --------
        >>> gr = pr.data.f1()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.max_disjoint(strand=False)
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        if strand is None:
            strand = self.strand_values_valid

        kwargs = {"strand": strand, "slack": slack}
        kwargs = fill_kwargs(kwargs)

        from pyranges.methods.max_disjoint import _max_disjoint

        dfs = pyrange_apply_single(_max_disjoint, self, **kwargs)

        return pr.PyRanges(dfs)

    def merge_overlaps(
        self,
        strand: Optional[bool] = None,
        count: bool = False,
        count_col: str = "Count",
        by: Optional[Union[List[str], str]] = None,
        slack: int = 0,
    ) -> "PyRanges":
        """Merge overlapping intervals into one.

        Parameters
        ----------
        strand : bool, default None, i.e. auto

            Only merge intervals on same strand.

        count : bool, default False

            Count intervals in each superinterval.

        count_col : str, default "Count"

            Name of column with counts.

        by : str or list of str, default None

            Only merge intervals with equal values in these columns.

        slack : int, default 0

            Allow this many nucleotides between each interval to merge.

        Returns
        -------
        PyRanges

            PyRanges with superintervals.

        Notes
        -----

        To avoid losing metadata, use cluster instead. If you want to perform a reduction function
        on the metadata, use pandas groupby.

        See Also
        --------

        PyRanges.cluster : annotate overlapping intervals with common ID

        Examples
        --------

        >>> gr = pr.data.ensembl_gtf()[["Feature", "gene_name"]]
        >>> gr
        +--------------+--------------+-----------+-----------+--------------+-------------+
        | Chromosome   | Feature      | Start     | End       | Strand       | gene_name   |
        | (category)   | (category)   | (int64)   | (int64)   | (category)   | (object)    |
        |--------------+--------------+-----------+-----------+--------------+-------------|
        | 1            | gene         | 11868     | 14409     | +            | DDX11L1     |
        | 1            | transcript   | 11868     | 14409     | +            | DDX11L1     |
        | 1            | exon         | 11868     | 12227     | +            | DDX11L1     |
        | 1            | exon         | 12612     | 12721     | +            | DDX11L1     |
        | ...          | ...          | ...       | ...       | ...          | ...         |
        | 1            | gene         | 1173055   | 1179555   | -            | TTLL10-AS1  |
        | 1            | transcript   | 1173055   | 1179555   | -            | TTLL10-AS1  |
        | 1            | exon         | 1179364   | 1179555   | -            | TTLL10-AS1  |
        | 1            | exon         | 1173055   | 1176396   | -            | TTLL10-AS1  |
        +--------------+--------------+-----------+-----------+--------------+-------------+
        Stranded PyRanges object has 2,446 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.merge_overlaps(count=True, count_col="Count")
        +--------------+-----------+-----------+--------------+-----------+
        | Chromosome   | Start     | End       | Strand       | Count     |
        | (category)   | (int64)   | (int64)   | (category)   | (int64)   |
        |--------------+-----------+-----------+--------------+-----------|
        | 1            | 11868     | 14409     | +            | 12        |
        | 1            | 29553     | 31109     | +            | 11        |
        | 1            | 52472     | 53312     | +            | 3         |
        | 1            | 57597     | 64116     | +            | 7         |
        | ...          | ...       | ...       | ...          | ...       |
        | 1            | 1062207   | 1063288   | -            | 4         |
        | 1            | 1070966   | 1074306   | -            | 10        |
        | 1            | 1081817   | 1116361   | -            | 319       |
        | 1            | 1173055   | 1179555   | -            | 4         |
        +--------------+-----------+-----------+--------------+-----------+
        Stranded PyRanges object has 62 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.merge_overlaps(by="Feature", count=True)
        +--------------+-----------+-----------+--------------+--------------+-----------+
        | Chromosome   | Start     | End       | Strand       | Feature      | Count     |
        | (category)   | (int64)   | (int64)   | (category)   | (category)   | (int64)   |
        |--------------+-----------+-----------+--------------+--------------+-----------|
        | 1            | 65564     | 65573     | +            | CDS          | 1         |
        | 1            | 69036     | 70005     | +            | CDS          | 2         |
        | 1            | 924431    | 924948    | +            | CDS          | 1         |
        | 1            | 925921    | 926013    | +            | CDS          | 11        |
        | ...          | ...       | ...       | ...          | ...          | ...       |
        | 1            | 1062207   | 1063288   | -            | transcript   | 1         |
        | 1            | 1070966   | 1074306   | -            | transcript   | 1         |
        | 1            | 1081817   | 1116361   | -            | transcript   | 19        |
        | 1            | 1173055   | 1179555   | -            | transcript   | 1         |
        +--------------+-----------+-----------+--------------+--------------+-----------+
        Stranded PyRanges object has 748 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.merge_overlaps(by=["Feature", "gene_name"], count=True)
        +--------------+-----------+-----------+--------------+--------------+-------------+-----------+
        | Chromosome   | Start     | End       | Strand       | Feature      | gene_name   | Count     |
        | (category)   | (int64)   | (int64)   | (category)   | (category)   | (object)    | (int64)   |
        |--------------+-----------+-----------+--------------+--------------+-------------+-----------|
        | 1            | 1020172   | 1020373   | +            | CDS          | AGRN        | 1         |
        | 1            | 1022200   | 1022462   | +            | CDS          | AGRN        | 2         |
        | 1            | 1034555   | 1034703   | +            | CDS          | AGRN        | 2         |
        | 1            | 1035276   | 1035324   | +            | CDS          | AGRN        | 4         |
        | ...          | ...       | ...       | ...          | ...          | ...         | ...       |
        | 1            | 347981    | 348366    | -            | transcript   | RPL23AP24   | 1         |
        | 1            | 1173055   | 1179555   | -            | transcript   | TTLL10-AS1  | 1         |
        | 1            | 14403     | 29570     | -            | transcript   | WASH7P      | 1         |
        | 1            | 185216    | 195411    | -            | transcript   | WASH9P      | 1         |
        +--------------+-----------+-----------+--------------+--------------+-------------+-----------+
        Stranded PyRanges object has 807 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        if strand is None:
            strand = self.strand_values_valid

        kwargs: Dict[str, Any] = {
            "strand": strand,
            "count": count,
            "by": by,
            "count_col": count_col,
            "slack": slack,
        }

        if not kwargs["by"]:
            kwargs["sparse"] = {"self": True}
            from pyranges.methods.merge import _merge

            df = pyrange_apply_single(_merge, self, **kwargs)
        else:
            kwargs["sparse"] = {"self": False}
            from pyranges.methods.merge import _merge_by

            df = pyrange_apply_single(_merge_by, self, **kwargs)

        return pr.PyRanges(df)

    def nearest(
        self,
        other: "PyRanges",
        strand_behavior: None = None,
        overlap: bool = True,
        how: Optional[str] = None,
        suffix: str = "_b",
        apply_strand_suffix: None = None,
    ) -> "PyRanges":
        """Find closest interval.

        Parameters
        ----------
        other : PyRanges

            PyRanges to find nearest interval in.

        strand_behavior : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        overlap : bool, default True

            Whether to include overlaps.

        how : {None, "upstream", "downstream"}, default None, i.e. both directions

            Whether to only look for nearest in one direction. Always with respect to the PyRanges
            it is called on.

        suffix : str, default "_b"

            Suffix to give columns with shared name in other.

        apply_strand_suffix : bool, default None

            If first pyranges is unstranded, but the second is not, the first will be given the strand column of the second.
            apply_strand_suffix makes the added strand column a regular data column instead by adding a suffix.

        Returns
        -------
        PyRanges

            A PyRanges with columns representing nearest interval horizontally appended.

        Notes
        -----

        A k_nearest also exists, but is less performant.

        See also
        --------

        PyRanges.new_position : give joined PyRanges new coordinates
        PyRanges.k_nearest : find k nearest intervals

        Examples
        --------

        >>> f1 = pr.PyRanges({'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 8, 5],
        ...                    'End': [6, 9, 7], 'Strand': ['+', '+', '-']})
        >>> f1
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         6 | +            |
        | chr1         |         8 |         9 | +            |
        | chr1         |         5 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f2 = pr.PyRanges({'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...                    'End': [2, 7], 'Strand': ['+', '-']})
        >>> f2
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         2 | +            |
        | chr1         |         6 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f1.nearest(f2)
        +--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------+
        | Chromosome   |     Start |       End | Strand       |   Start_b |     End_b | Strand_b     |   Distance |
        | (category)   |   (int64) |   (int64) | (category)   |   (int64) |   (int64) | (category)   |    (int64) |
        |--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------|
        | chr1         |         3 |         6 | +            |         6 |         7 | -            |          1 |
        | chr1         |         8 |         9 | +            |         6 |         7 | -            |          2 |
        | chr1         |         5 |         7 | -            |         6 |         7 | -            |          0 |
        +--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 3 rows and 8 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> f1.nearest(f2, how="upstream")
        +--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------+
        | Chromosome   |     Start |       End | Strand       |   Start_b |     End_b | Strand_b     |   Distance |
        | (category)   |   (int64) |   (int64) | (category)   |   (int64) |   (int64) | (category)   |    (int64) |
        |--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------|
        | chr1         |         3 |         6 | +            |         1 |         2 | +            |          2 |
        | chr1         |         8 |         9 | +            |         6 |         7 | -            |          2 |
        | chr1         |         5 |         7 | -            |         6 |         7 | -            |          0 |
        +--------------+-----------+-----------+--------------+-----------+-----------+--------------+------------+
        Stranded PyRanges object has 3 rows and 8 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        from pyranges.methods.nearest import _nearest

        kwargs = {
            "strand_behavior": strand_behavior,
            "how": how,
            "overlap": overlap,
            "suffix": suffix,
            "apply_strand_suffix": apply_strand_suffix,
        }
        kwargs = fill_kwargs(kwargs)
        if kwargs.get("how") in "upstream downstream".split():
            assert other.strand_values_valid, "If doing upstream or downstream nearest, other pyranges must be stranded"

        dfs = pyrange_apply(_nearest, self, other, **kwargs)
        gr = pr.PyRanges(dfs)

        if not self.strand_values_valid and other.strand_values_valid:
            if apply_strand_suffix is None:
                import sys

                print(
                    "join: Strand data from other will be added as strand data to self.\nIf this is undesired use the flag apply_strand_suffix=False.\nTo turn off the warning set apply_strand_suffix to True or False.",
                    file=sys.stderr,
                )
            elif apply_strand_suffix:
                gr.columns = gr.columns.str.replace(
                    "Strand", "Strand" + kwargs["suffix"]
                )

        return gr

    def new_position(
        self, new_pos: str, columns: Optional[Tuple[str, str, str, str]] = None
    ) -> "PyRanges":
        """Give new position.

        The operation join produces a PyRanges with two pairs of start coordinates and two pairs of
        end coordinates. This operation uses these to give the PyRanges a new position.

        Parameters
        ----------
        new_pos : {"union", "intersection", "swap"}

           Change of coordinates.

        columns : Optional[tuple of str], default None, i.e. auto

           The name of the coordinate columns. By default, uses the two first columns containing
           "Start" and the two first columns containing "End".

        See Also
        --------

        PyRanges.join : combine two PyRanges horizontally with SQL-style joins.

        Returns
        -------
        PyRanges

            PyRanges with new coordinates.

        Examples
        --------

        >>> gr = pr.PyRanges({'Chromosome': ['chr1', 'chr1', 'chr1'],
        ...                    'Start': [3, 8, 5], 'End': [6, 9, 7]})
        >>> gr
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int64) |   (int64) |
        |--------------+-----------+-----------|
        | chr1         |         3 |         6 |
        | chr1         |         8 |         9 |
        | chr1         |         5 |         7 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr2 = pr.PyRanges({'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...                     'End': [4, 7]})
        >>> gr2
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int64) |   (int64) |
        |--------------+-----------+-----------|
        | chr1         |         1 |         4 |
        | chr1         |         6 |         7 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> j = gr.join_overlaps(gr2)
        >>> j
        +--------------+-----------+-----------+-----------+-----------+
        | Chromosome   |     Start |       End |   Start_b |     End_b |
        | (category)   |   (int64) |   (int64) |   (int64) |   (int64) |
        |--------------+-----------+-----------+-----------+-----------|
        | chr1         |         3 |         6 |         1 |         4 |
        | chr1         |         5 |         7 |         6 |         7 |
        +--------------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> j.new_position("swap")
        +--------------+-----------+-----------+-----------+-----------+
        | Chromosome   |     Start |       End |   Start_b |     End_b |
        | (category)   |   (int64) |   (int64) |   (int64) |   (int64) |
        |--------------+-----------+-----------+-----------+-----------|
        | chr1         |         1 |         4 |         3 |         6 |
        | chr1         |         6 |         7 |         5 |         7 |
        +--------------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> j.new_position("union").mp()
        +--------------------+-----------+-----------+
        | - Position -       |   Start_b |     End_b |
        | (Multiple types)   |   (int64) |   (int64) |
        |--------------------+-----------+-----------|
        | chr1 1-6           |         1 |         4 |
        | chr1 5-7           |         6 |         7 |
        +--------------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> j.new_position("intersection").mp()
        +--------------------+-----------+-----------+
        | - Position -       |   Start_b |     End_b |
        | (Multiple types)   |   (int64) |   (int64) |
        |--------------------+-----------+-----------|
        | chr1 1-4           |         1 |         4 |
        | chr1 6-7           |         6 |         7 |
        +--------------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> j2 = pr.PyRanges({"Chromosome": [1], "Start": [3],
        ...                   "End": [4], "A": [1], "B": [3], "C": [2], "D": [5]})
        >>> j2
        +--------------+-----------+-----------+-----------+-----------+-----------+-----------+
        |   Chromosome |     Start |       End |         A |         B |         C |         D |
        |   (category) |   (int64) |   (int64) |   (int64) |   (int64) |   (int64) |   (int64) |
        |--------------+-----------+-----------+-----------+-----------+-----------+-----------|
        |            1 |         3 |         4 |         1 |         3 |         2 |         5 |
        +--------------+-----------+-----------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 1 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> j2.new_position("intersection", ("A", "B", "C", "D"))
        +--------------+-----------+-----------+-----------+-----------+-----------+-----------+
        |   Chromosome |     Start |       End |         A |         B |         C |         D |
        |   (category) |   (int64) |   (int64) |   (int64) |   (int64) |   (int64) |   (int64) |
        |--------------+-----------+-----------+-----------+-----------+-----------+-----------|
        |            1 |         2 |         3 |         1 |         3 |         2 |         5 |
        +--------------+-----------+-----------+-----------+-----------+-----------+-----------+
        Unstranded PyRanges object has 1 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        """

        from pyranges.methods.new_position import _new_position

        if self.empty:
            return self

        kwargs: Dict[str, Any] = {
            "strand": None,
            "sparse": {"self": False},
            "new_pos": new_pos,
        }

        if columns is None:
            start1, start2 = self.columns[self.columns.str.contains("Start")][:2]
            end1, end2 = self.columns[self.columns.str.contains("End")][:2]
            columns = (start1, end1, start2, end2)

        kwargs["columns"] = columns

        kwargs = fill_kwargs(kwargs)

        dfs = pyrange_apply_single(_new_position, self, **kwargs)

        return pr.PyRanges(dfs)

    def overlap(
        self,
        other: "PyRanges",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        how: VALID_OVERLAP_TYPE = "first",
        invert: bool = False,
    ) -> "PyRanges":
        """Return overlapping intervals.

        Returns the intervals in self which overlap with those in other.

        Parameters
        ----------
        other : PyRanges

            PyRanges to find overlaps with.

        strand_behavior : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        how : {"first", "containment", False, None}, default "first"

            What intervals to report. By default, reports every interval in self with overlap once.
            "containment" reports all intervals where the overlapping is contained within it.

        invert : bool, default False

            Whether to return the intervals without overlaps.

        Returns
        -------
        PyRanges

            A PyRanges with overlapping intervals.

        See also
        --------

        PyRanges.intersect : report overlapping subintervals
        PyRanges.set_intersect : set-intersect PyRanges

        Examples
        --------

        >>> gr = pr.PyRanges({"Chromosome": ["chr1", "chr1", "chr2", "chr1", "chr3"], "Start": [1, 1, 4, 10, 0],
        ...                    "End": [3, 3, 9, 11, 1], "ID": ["A", "a", "b", "c", "d"]})
        >>> gr2 = pr.PyRanges({"Chromosome": ["chr1", "chr1", "chr2"], "Start": [2, 2, 1], "End": [3, 9, 10]})
        >>> gr.overlap(gr2, how="first")
        Chromosome      Start      End  ID
        object          int64    int64  object
        ------------  -------  -------  --------
        chr1                1        3  A
        chr1                1        3  a
        chr2                4        9  b

        >>> gr.overlap(gr2, how="containment")
        Chromosome      Start      End  ID
        object          int64    int64  object
        ------------  -------  -------  --------
        chr2                4        9  b

        >>> gr.overlap(gr2, invert=True)
        Chromosome      Start      End  ID
        object          int64    int64  object
        ------------  -------  -------  --------
        chr1               10       11  c
        chr3                0        1  d
        """
        return PyRanges(
            pyrange_apply(
                _overlap,
                self,
                other,
                invert=invert,
                how=how,
                strand_behavior=strand_behavior,
            )
        )

    def sample(self, *args, **kwargs) -> "PyRanges":
        return PyRanges(super().sample(*args, **kwargs))

    def set_intersect(
        self,
        other: "PyRanges",
        strand_behavior: Optional[str] = None,
        how: Optional[str] = None,
        new_pos: bool = False,
    ) -> "PyRanges":
        """Return set-theoretical intersection.

        Like intersect, but both PyRanges are merged first.

        Parameters
        ----------
        other : PyRanges

            PyRanges to set-intersect.

        strand_behavior : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        how : {None, "first", "last", "containment"}, default None, i.e. all

            What intervals to report. By default, reports all overlapping intervals. "containment"
            reports intervals where the overlapping is contained within it.

        Returns
        -------
        PyRanges

            A PyRanges with overlapping subintervals.

        See also
        --------

        PyRanges.intersect : find overlapping subintervals
        PyRanges.overlap : report overlapping intervals

        Examples
        --------

        >>> gr = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [1, 4, 10],
        ...                    "End": [3, 9, 11], "ID": ["a", "b", "c"]})
        >>> gr
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | ID         |
        | (category)   |   (int64) |   (int64) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         1 |         3 | a          |
        | chr1         |         4 |         9 | b          |
        | chr1         |        10 |        11 | c          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr2 = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [2, 2, 9], "End": [3, 9, 10]})
        >>> gr2
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int64) |   (int64) |
        |--------------+-----------+-----------|
        | chr1         |         2 |         3 |
        | chr1         |         2 |         9 |
        | chr1         |         9 |        10 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.set_intersect(gr2)
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int64) |   (int64) |
        |--------------+-----------+-----------|
        | chr1         |         2 |         3 |
        | chr1         |         4 |         9 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        In this simple unstranded case, this is the same as the below:

        >>> gr.merge_overlaps().intersect(gr2.merge_overlaps())
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int64) |   (int64) |
        |--------------+-----------+-----------|
        | chr1         |         2 |         3 |
        | chr1         |         4 |         9 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.set_intersect(gr2, how="containment")
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int64) |   (int64) |
        |--------------+-----------+-----------|
        | chr1         |         4 |         9 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 1 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        kwargs = {
            "strand_behavior": strand_behavior,
            "how": how,
            "new_pos": new_pos,
        }
        kwargs = fill_kwargs(kwargs)
        strand = True if strand_behavior else False
        self_clusters = self.merge_overlaps(strand=strand)
        other_clusters = other.merge_overlaps(strand=strand)
        dfs = pyrange_apply(_intersection, self_clusters, other_clusters, **kwargs)

        return pr.PyRanges(dfs)

    def set_union(self, other: "PyRanges", strand_behavior: None = None) -> "PyRanges":
        """Return set-theoretical union.

        Parameters
        ----------
        other : PyRanges

            PyRanges to do union with.

        strand_behavior : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        Returns
        -------
        PyRanges

            A PyRanges with the union of intervals.

        See also
        --------

        PyRanges.set_intersect : set-theoretical intersection
        PyRanges.overlap : report overlapping intervals

        Examples
        --------

        >>> gr = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [1, 4, 10],
        ...                    "End": [3, 9, 11], "ID": ["a", "b", "c"]})
        >>> gr
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | ID         |
        | (category)   |   (int64) |   (int64) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         1 |         3 | a          |
        | chr1         |         4 |         9 | b          |
        | chr1         |        10 |        11 | c          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr2 = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [2, 2, 9], "End": [3, 9, 10]})
        >>> gr2
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int64) |   (int64) |
        |--------------+-----------+-----------|
        | chr1         |         2 |         3 |
        | chr1         |         2 |         9 |
        | chr1         |         9 |        10 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.set_union(gr2)
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int64) |   (int64) |
        |--------------+-----------+-----------|
        | chr1         |         1 |        11 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 1 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        if self.empty and other.empty:
            return pr.PyRanges()

        strand = True if strand_behavior else False

        if not strand:
            self = self.remove_strand()
            other = other.remove_strand()

        if strand_behavior == "opposite" and len(other):
            other = other.copy()
            other.Strand = other.Strand.replace({"+": "-", "-": "+"})

        gr = pr.concat([self, other], strand)

        gr = gr.merge_overlaps(strand=strand)

        return gr

    def sort(self, by: Optional[str] = None) -> "PyRanges":
        """Sort by position or columns.

        Parameters
        ----------
        by : str or list of str, default None

            Column(s) to sort by. Default is Start and End.
            Special value "5" can be provided to sort by 5': intervals on + strand are sorted in ascending order, while
            those on - strand are sorted in descending order.

        Note
        ----

        Since a PyRanges contains multiple pd.DataFrames, the sorting only happens within dataframes.

        Returns
        -------
        PyRanges

            Sorted PyRanges

        See Also
        --------

        pyranges.multioverlap : find overlaps with multiple PyRanges

        Examples
        --------

        >>> p  = pr.PyRanges({"Chromosome": [1, 1, 1, 1, 1, 1],
        ...                    "Strand": ["+", "+", "-", "-", "+", "+"],
        ...                    "Start": [40, 1, 10, 70, 140, 160],
        ...                    "End": [60, 11, 25, 80, 152, 190],
        ...                    "transcript_id":["t3", "t3", "t2", "t2", "t1", "t1"] })

        By default, intervals are sorted by position:

        >>> p.sort()
        +--------------+--------------+-----------+-----------+-----------------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |
        |--------------+--------------+-----------+-----------+-----------------|
        |            1 | +            |         1 |        11 | t3              |
        |            1 | +            |        40 |        60 | t3              |
        |            1 | +            |       140 |       152 | t1              |
        |            1 | +            |       160 |       190 | t1              |
        |            1 | -            |        10 |        25 | t2              |
        |            1 | -            |        70 |        80 | t2              |
        +--------------+--------------+-----------+-----------+-----------------+
        Stranded PyRanges object has 6 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        (Note how sorting takes place within Chromosome-Strand pairs.)

        To sort according to a specified column:

        >>> p.sort(by='transcript_id')
        +--------------+--------------+-----------+-----------+-----------------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |
        |--------------+--------------+-----------+-----------+-----------------|
        |            1 | +            |       140 |       152 | t1              |
        |            1 | +            |       160 |       190 | t1              |
        |            1 | +            |        40 |        60 | t3              |
        |            1 | +            |         1 |        11 | t3              |
        |            1 | -            |        10 |        25 | t2              |
        |            1 | -            |        70 |        80 | t2              |
        +--------------+--------------+-----------+-----------+-----------------+
        Stranded PyRanges object has 6 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        If the special value "5" is provided, intervals are sorted
        according to their five-prime end:

        >>> p.sort("5")
        +--------------+--------------+-----------+-----------+-----------------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |
        |--------------+--------------+-----------+-----------+-----------------|
        |            1 | +            |         1 |        11 | t3              |
        |            1 | +            |        40 |        60 | t3              |
        |            1 | +            |       140 |       152 | t1              |
        |            1 | +            |       160 |       190 | t1              |
        |            1 | -            |        70 |        80 | t2              |
        |            1 | -            |        10 |        25 | t2              |
        +--------------+--------------+-----------+-----------+-----------------+
        Stranded PyRanges object has 6 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        """

        from pyranges.methods.sort import _sort

        kwargs = {"strand": self.strand_values_valid, "sparse": {"self": False}}
        if by:
            assert "5" not in by or (
                ((type(by) is str and by == "5") or (type(by) is not str and "5" in by))
                and self.strand_values_valid
            ), "Only stranded PyRanges can be sorted by 5'! "
            kwargs["by"] = by

        kwargs = fill_kwargs(kwargs)
        return pr.PyRanges(pyrange_apply_single(_sort, self, **kwargs))

    def sp(self, n=30, formatting=None):
        """Sort on location and print.

        See Also
        --------

        PyRanges.print : print PyRanges."""

        print(tostring(self, n=n, sort=True, formatting=formatting))

    def spc(self, n=30, formatting=None):
        """Sort on location, print and return self.

        See Also
        --------

        PyRanges.print : print PyRanges."""

        print(tostring(self, n=n, sort=True, formatting=formatting))

        return self

    def slack(self, slack):
        """Deprecated: this function has been moved to Pyranges.extend"""
        return self.extend(slack)

    def spliced_subsequence(
        self,
        start: int = 0,
        end: Optional[int] = None,
        by: Optional[str] = None,
        strand: Optional[bool] = None,
        **kwargs,
    ) -> "PyRanges":
        """Get subsequences of the intervals, using coordinates mapping to spliced transcripts (without introns)

        The returned intervals are subregions of self, cut according to specifications.
        Start and end are relative to the 5' end: 0 means the leftmost nucleotide for + strand
        intervals, while it means the rightmost one for - strand.
        This method also allows to manipulate groups of intervals (e.g. exons belonging to same transcripts)
        through the 'by' argument. When using it, start and end refer to the spliced transcript coordinates,
        meaning that introns are ignored in the count.

        Parameters
        ----------

        start : int
            Start of subregion, 0-based and included, counting from the 5' end.
            Use a negative int to count from the 3'  (e.g. -1 is the last nucleotide)

        end : int, default None
            End of subregion, 0-based and excluded, counting from the 5' end.
            Use a negative int to count from the 3'  (e.g. -1 is the last nucleotide)
            If None, the existing 3' end is returned.

        by : list of str, default None
            intervals are grouped by this/these ID column(s) beforehand, e.g. exons belonging to same transcripts


        strand : bool, default None, i.e. auto
            Whether strand is considered when interpreting the start and end arguments of this function.
            If True, counting is from the 5' end, which is the leftmost coordinate for + strand and the rightmost for - strand.
            If False, all intervals are processed like they reside on the + strand.
            If None (default), strand is considered if the PyRanges is stranded.

        Returns
        -------

        PyRanges
            Subregion of self, subsequenced as specified by arguments

        Note
        ----
        If the request goes out of bounds (e.g. requesting 100 nts for a 90nt region), only the existing portion is returned

        See also
        --------
        subsequence : analogous to this method, but input coordinates refer to the unspliced transcript

        Examples
        --------
        >>> p  = pr.PyRanges({"Chromosome": [1, 1, 2, 2, 3],
        ...                   "Strand": ["+", "+", "-", "-", "+"],
        ...                   "Start": [1, 40, 10, 70, 140],
        ...                   "End": [11, 60, 25, 80, 152],
        ...                   "transcript_id":["t1", "t1", "t2", "t2", "t3"] })
        >>> p
        +--------------+--------------+-----------+-----------+-----------------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |
        |--------------+--------------+-----------+-----------+-----------------|
        |            1 | +            |         1 |        11 | t1              |
        |            1 | +            |        40 |        60 | t1              |
        |            2 | -            |        10 |        25 | t2              |
        |            2 | -            |        70 |        80 | t2              |
        |            3 | +            |       140 |       152 | t3              |
        +--------------+--------------+-----------+-----------+-----------------+
        Stranded PyRanges object has 5 rows and 5 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

          # Get the first 15 nucleotides of *each spliced transcript*, grouping exons by transcript_id:
        >>> p.spliced_subsequence(0, 15, by='transcript_id')
        +--------------+--------------+-----------+-----------+-----------------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |
        |--------------+--------------+-----------+-----------+-----------------|
        |            1 | +            |         1 |        11 | t1              |
        |            1 | +            |        40 |        45 | t1              |
        |            2 | -            |        70 |        80 | t2              |
        |            2 | -            |        20 |        25 | t2              |
        |            3 | +            |       140 |       152 | t3              |
        +--------------+--------------+-----------+-----------+-----------------+
        Stranded PyRanges object has 5 rows and 5 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

          # Get the last 20 nucleotides of each spliced transcript:
        >>> p.spliced_subsequence(-20, by='transcript_id')
        +--------------+--------------+-----------+-----------+-----------------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |
        |--------------+--------------+-----------+-----------+-----------------|
        |            1 | +            |        40 |        60 | t1              |
        |            2 | -            |        70 |        75 | t2              |
        |            2 | -            |        10 |        25 | t2              |
        |            3 | +            |       140 |       152 | t3              |
        +--------------+--------------+-----------+-----------+-----------------+
        Stranded PyRanges object has 4 rows and 5 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

          # Get region from 25 to 60 of each spliced transcript, or their existing subportion:
        >>> p.spliced_subsequence(25, 60, by='transcript_id')
        +--------------+--------------+-----------+-----------+-----------------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |
        |--------------+--------------+-----------+-----------+-----------------|
        |            1 | +            |        55 |        60 | t1              |
        +--------------+--------------+-----------+-----------+-----------------+
        Stranded PyRanges object has 1 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

          # Get region of each spliced transcript which excludes their first and last 3 nucleotides:
        >>> p.spliced_subsequence(3, -3, by='transcript_id')
        +--------------+--------------+-----------+-----------+-----------------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |
        |--------------+--------------+-----------+-----------+-----------------|
        |            1 | +            |         4 |        11 | t1              |
        |            1 | +            |        40 |        57 | t1              |
        |            2 | -            |        70 |        77 | t2              |
        |            2 | -            |        13 |        25 | t2              |
        |            3 | +            |       143 |       149 | t3              |
        +--------------+--------------+-----------+-----------+-----------------+
        Stranded PyRanges object has 5 rows and 5 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        from pyranges.methods.spliced_subsequence import _spliced_subseq

        if strand and not self.strand_values_valid:
            raise Exception(
                "spliced_subsequence: you can use strand=True only for stranded PyRanges!"
            )

        if strand is None:
            strand = True if self.strand_values_valid else False

        kwargs.update({"strand": strand, "by": by, "start": start, "end": end})
        kwargs = fill_kwargs(kwargs)

        if not strand:
            sorted_p = self.sort()
        else:
            sorted_p = self.sort("5")

        result = pyrange_apply_single(_spliced_subseq, sorted_p, **kwargs)

        return pr.PyRanges(result)

    def split(self, strand: Optional[bool] = None, between: bool = False) -> "PyRanges":
        """Split into non-overlapping intervals.

        Parameters
        ----------
        strand : Optional[bool], default None, i.e. auto

            Whether to ignore strand information if PyRanges is stranded.

        between : bool, default False

            Include lengths between intervals.

        Returns
        -------
        PyRanges

            PyRanges with intervals split at overlap points.

        See Also
        --------

        pyranges.multioverlap : find overlaps with multiple PyRanges

        Examples
        --------

        >>> d = {'Chromosome': ['chr1', 'chr1', 'chr1', 'chr1'], 'Start': [3, 5, 5, 11],
        ...       'End': [6, 9, 7, 12], 'Strand': ['+', '+', '-', '-']}
        >>> gr = pr.PyRanges(d)
        >>> gr
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         6 | +            |
        | chr1         |         5 |         9 | +            |
        | chr1         |         5 |         7 | -            |
        | chr1         |        11 |        12 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 4 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.split()
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         5 | +            |
        | chr1         |         5 |         6 | +            |
        | chr1         |         6 |         9 | +            |
        | chr1         |         5 |         7 | -            |
        | chr1         |        11 |        12 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 5 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.split(between=True)
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         3 |         5 | +            |
        | chr1         |         5 |         6 | +            |
        | chr1         |         6 |         9 | +            |
        | chr1         |         5 |         7 | -            |
        | chr1         |         7 |        11 | -            |
        | chr1         |        11 |        12 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 6 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.split(strand=False)
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int64) |   (int64) |
        |--------------+-----------+-----------|
        | chr1         |         3 |         5 |
        | chr1         |         5 |         6 |
        | chr1         |         6 |         7 |
        | chr1         |         7 |         9 |
        | chr1         |        11 |        12 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 5 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.split(strand=False, between=True)
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int64) |   (int64) |
        |--------------+-----------+-----------|
        | chr1         |         3 |         5 |
        | chr1         |         5 |         6 |
        | chr1         |         6 |         7 |
        | chr1         |         7 |         9 |
        | chr1         |         9 |        11 |
        | chr1         |        11 |        12 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 6 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        if strand is None:
            strand = self.strand_values_valid

        kwargs = fill_kwargs({"strand": strand})

        from pyranges.methods.split import _split

        df = pyrange_apply_single(_split, self, **kwargs)

        split = pr.PyRanges(df)
        if not between:
            strand_behavior: Union[str, bool] = "same" if strand else False
            split = split.overlap(self, strand_behavior=strand_behavior)

        return split

    @property
    def strand_values_valid(self) -> bool:
        """Whether PyRanges has (valid) strand info.

        Note
        ----

        A PyRanges can have invalid values in the Strand-column. It is not considered stranded.

        See Also
        --------

        PyRanges.strands : return the strands

        Examples
        --------

        >>> d =  {'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...       'End': [5, 8], 'Strand': ['+', '.']}
        >>> gr = pr.PyRanges(d)
        >>> gr
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         5 | +            |
        | chr1         |         6 |         8 | .            |
        +--------------+-----------+-----------+--------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        Considered unstranded due to these Strand values: '.'

        >>> gr.strand_values_valid  # invalid strand value: '.'
        False

        >>> "Strand" in gr.columns
        True
        """
        if STRAND_COL not in self.columns and len(self) > 0:
            return False
        else:
            return self[STRAND_COL].isin(VALID_GENOMIC_STRAND_INFO).all()

    @property
    def strands(self) -> List[Union[Any, str]]:
        """Return strands.

        Notes
        -----

        If the strand-column contains an invalid value, [] is returned.

        See Also
        --------

        PyRanges.stranded : whether has valid strand info

        Examples
        --------
        >>> d =  {'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...       'End': [5, 8], 'Strand': ['+', '.']}
        >>> gr = pr.PyRanges(d)
        >>> gr
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         5 | +            |
        | chr1         |         6 |         8 | .            |
        +--------------+-----------+-----------+--------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        Considered unstranded due to these Strand values: '.'

        >>> gr.strands
        []

        >>> gr.Strand.drop_duplicates().to_list()
        ['+', '.']

        >>> gr.Strand = ["+", "-"]
        >>> gr.strands
        ['+', '-']
        """

        if not self.strand_values_valid:
            return []

        return natsorted(set([k[1] for k in self.keys()]))

    def subset(
        self, f: Callable, strand: Optional[bool] = None, **kwargs
    ) -> "PyRanges":
        """Return a subset of the rows.

        Parameters
        ----------
        f : function
            Function which returns boolean pd.Series equal to length of df.

        strand : bool, default None, i.e. auto

            Whether to do operations on chromosome/strand pairs or chromosomes. If None, will use
            chromosome/strand pairs if the PyRanges is stranded.

        **kwargs
            Additional keyword arguments to pass as keyword arguments to `f`

        Notes
        -----

        PyRanges can also be subsetted directly with a boolean pd.Series. This function is slightly
        faster, but more cumbersome.

        Returns
        -------
        PyRanges

            PyRanges subset on rows.

        Examples
        --------

        >>> gr = pr.data.f1()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         3 |         6 | interval1  |         0 | +            |
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.subset(lambda df: df.Start > 4)
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        Also possible:

        >>> gr[gr.Start > 4]
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         |         8 |         9 | interval3  |         0 | +            |
        | chr1         |         5 |         7 | interval2  |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        kwargs = fill_kwargs(kwargs)

        if strand is None:
            strand = self.strand_values_valid

        if self.strand_values_valid and not strand:
            self = self.remove_strand()

        kwargs.update({"strand": strand})

        result = pyrange_apply_single(f, self, **kwargs)

        if not result:
            return pr.PyRanges()

        first_result = next(iter(result.values()))

        assert (
            first_result.dtype == bool
        ), "result of subset function must be bool, but is {}".format(
            first_result.dtype
        )

        return self[result]

    def subsequence(
        self,
        start: int = 0,
        end: Optional[int] = None,
        by: Optional[str] = None,
        strand: Optional[bool] = None,
        **kwargs,
    ) -> "PyRanges":
        """Get subsequences of the intervals.

        The returned intervals are subregions of self, cut according to specifications.
        Start and end are relative to the 5' end: 0 means the leftmost nucleotide for + strand
        intervals, while it means the rightmost one for - strand.
        This method also allows to manipulate groups of intervals (e.g. exons belonging to same transcripts)
        through the 'by' argument. When using it, start and end refer to the unspliced transcript coordinates,
        meaning that introns are included in the count.

        Parameters
        ----------

        start : int
            Start of subregion, 0-based and included, counting from the 5' end.
            Use a negative int to count from the 3'  (e.g. -1 is the last nucleotide)

        end : int, default None
            End of subregion, 0-based and excluded, counting from the 5' end.
            Use a negative int to count from the 3'  (e.g. -1 is the last nucleotide)

            If None, the existing 3' end is returned.

        by : list of str, default None
            intervals are grouped by this/these ID column(s) beforehand, e.g. exons belonging to same transcripts

        strand : bool, default None, i.e. auto
            Whether strand is considered when interpreting the start and end arguments of this function.
            If True, counting is from the 5' end, which is the leftmost coordinate for + strand and the rightmost for - strand.
            If False, all intervals are processed like they reside on the + strand.
            If None (default), strand is considered if the PyRanges is stranded.

        Returns
        -------

        PyRanges
            Subregion of self, subsequenced as specified by arguments

        Note
        ----
        If the request goes out of bounds (e.g. requesting 100 nts for a 90nt region), only the existing portion is returned

        See also
        --------
        spliced_subsequence : analogous to this method, but intronic regions are not counted, so that input coordinates refer to the spliced transcript


        Examples
        --------
        >>> p  = pr.PyRanges({"Chromosome": [1, 1, 2, 2, 3],
        ...                   "Strand": ["+", "+", "-", "-", "+"],
        ...                   "Start": [1, 40, 2, 30, 140],
        ...                   "End": [20, 60, 13, 45, 155],
        ...                   "transcript_id":["t1", "t1", "t2", "t2", "t3"] })
        >>> p
        +--------------+--------------+-----------+-----------+-----------------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |
        |--------------+--------------+-----------+-----------+-----------------|
        |            1 | +            |         1 |        20 | t1              |
        |            1 | +            |        40 |        60 | t1              |
        |            2 | -            |         2 |        13 | t2              |
        |            2 | -            |        30 |        45 | t2              |
        |            3 | +            |       140 |       155 | t3              |
        +--------------+--------------+-----------+-----------+-----------------+
        Stranded PyRanges object has 5 rows and 5 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        # Get the first 10 nucleotides (at the 5') of *each interval* (each line of the dataframe):
        >>> p.subsequence(0, 10)
        +--------------+--------------+-----------+-----------+-----------------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |
        |--------------+--------------+-----------+-----------+-----------------|
        |            1 | +            |         1 |        11 | t1              |
        |            1 | +            |        40 |        50 | t1              |
        |            2 | -            |         3 |        13 | t2              |
        |            2 | -            |        35 |        45 | t2              |
        |            3 | +            |       140 |       150 | t3              |
        +--------------+--------------+-----------+-----------+-----------------+
        Stranded PyRanges object has 5 rows and 5 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        # Get the first 10 nucleotides of *each transcript*, grouping exons by transcript_id:
        >>> p.subsequence(0, 10, by='transcript_id')
        +--------------+--------------+-----------+-----------+-----------------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |
        |--------------+--------------+-----------+-----------+-----------------|
        |            1 | +            |         1 |        11 | t1              |
        |            2 | -            |        35 |        45 | t2              |
        |            3 | +            |       140 |       150 | t3              |
        +--------------+--------------+-----------+-----------+-----------------+
        Stranded PyRanges object has 3 rows and 5 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        # Get the last 20 nucleotides of each transcript:
        >>> p.subsequence(-20, by='transcript_id')
        +--------------+--------------+-----------+-----------+-----------------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |
        |--------------+--------------+-----------+-----------+-----------------|
        |            1 | +            |        40 |        60 | t1              |
        |            2 | -            |         2 |        13 | t2              |
        |            3 | +            |       140 |       155 | t3              |
        +--------------+--------------+-----------+-----------+-----------------+
        Stranded PyRanges object has 3 rows and 5 columns from 3 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        # Get region from 30 to 330 of each transcript, or their existing subportion:
        >>> p.subsequence(30, 300, by='transcript_id')
        +--------------+--------------+-----------+-----------+-----------------+
        |   Chromosome | Strand       |     Start |       End | transcript_id   |
        |   (category) | (category)   |   (int64) |   (int64) | (object)        |
        |--------------+--------------+-----------+-----------+-----------------|
        |            1 | +            |        40 |        60 | t1              |
        |            2 | -            |         2 |        13 | t2              |
        +--------------+--------------+-----------+-----------+-----------------+
        Stranded PyRanges object has 2 rows and 5 columns from 2 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """
        from pyranges.methods.subsequence import _subseq

        if strand is None:
            strand = True if self.strand_values_valid else False

        kwargs.update({"strand": strand, "by": by, "start": start, "end": end})
        kwargs = fill_kwargs(kwargs)

        result = pyrange_apply_single(_subseq, self, **kwargs)

        return pr.PyRanges(result)

    def range_subtract(
        self, other: "PyRanges", strand_behavior: None = None
    ) -> "PyRanges":
        """Subtract intervals.

        Parameters
        ----------
        strand_behavior : {None, "same", "opposite", False}, default None, i.e. auto

            Whether to compare PyRanges on the same strand, the opposite or ignore strand
            information. The default, None, means use "same" if both PyRanges are strande,
            otherwise ignore the strand information.

        See Also
        --------
        pyranges.PyRanges.overlap : use with invert=True to return all intervals without overlap

        Examples
        --------

        >>> gr = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [1, 4, 10],
        ...                    "End": [3, 9, 11], "ID": ["a", "b", "c"]})
        >>> gr2 = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [2, 2, 9], "End": [3, 9, 10]})
        >>> gr
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | ID         |
        | (category)   |   (int64) |   (int64) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         1 |         3 | a          |
        | chr1         |         4 |         9 | b          |
        | chr1         |        10 |        11 | c          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 3 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr2
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int64) |   (int64) |
        |--------------+-----------+-----------|
        | chr1         |         2 |         3 |
        | chr1         |         2 |         9 |
        | chr1         |         9 |        10 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 3 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.range_subtract(gr2)
        +--------------+-----------+-----------+------------+
        | Chromosome   |     Start |       End | ID         |
        | (category)   |   (int64) |   (int64) | (object)   |
        |--------------+-----------+-----------+------------|
        | chr1         |         1 |         2 | a          |
        | chr1         |        10 |        11 | c          |
        +--------------+-----------+-----------+------------+
        Unstranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        from pyranges.methods.subtraction import _subtraction

        kwargs = {
            "strand_behavior": strand_behavior,
            "sparse": {"self": False, "other": True},
        }
        kwargs = fill_kwargs(kwargs)

        strand = True if strand_behavior else False
        other_clusters = other.merge_overlaps(strand=strand)

        _self = self.copy()

        _self = _self.count_overlaps(
            other_clusters, strand_behavior=strand_behavior, overlap_col="__num__"
        )

        result = pyrange_apply(_subtraction, _self, other_clusters, **kwargs)

        return PyRanges(pr.PyRanges(result).drop("__num__"))

    def summary(
        self, to_stdout: bool = True, return_df: bool = False
    ) -> Optional[pd.DataFrame]:
        """Return info.

        Count refers to the number of intervals, the rest to the lengths.

        The column "pyrange" describes the data as is. "coverage_forward" and "coverage_reverse"
        describe the data after strand-specific merging of overlapping intervals.
        "coverage_unstranded" describes the data after merging, without considering the strands.

        The row "count" is the number of intervals and "sum" is their total length. The rest describe the lengths of the
        intervals.

        Parameters
        ----------
        to_stdout : bool, default True

            Print summary.

        return_df : bool, default False

            Return df with summary.

        Returns
        -------
            None or pd.DataFrame with summary.


        Examples
        --------

        >>> gr = pr.data.ensembl_gtf()[["Feature", "gene_id"]]
        >>> gr
        +--------------+--------------+-----------+-----------+--------------+-----------------+
        | Chromosome   | Feature      | Start     | End       | Strand       | gene_id         |
        | (category)   | (category)   | (int64)   | (int64)   | (category)   | (object)        |
        |--------------+--------------+-----------+-----------+--------------+-----------------|
        | 1            | gene         | 11868     | 14409     | +            | ENSG00000223972 |
        | 1            | transcript   | 11868     | 14409     | +            | ENSG00000223972 |
        | 1            | exon         | 11868     | 12227     | +            | ENSG00000223972 |
        | 1            | exon         | 12612     | 12721     | +            | ENSG00000223972 |
        | ...          | ...          | ...       | ...       | ...          | ...             |
        | 1            | gene         | 1173055   | 1179555   | -            | ENSG00000205231 |
        | 1            | transcript   | 1173055   | 1179555   | -            | ENSG00000205231 |
        | 1            | exon         | 1179364   | 1179555   | -            | ENSG00000205231 |
        | 1            | exon         | 1173055   | 1176396   | -            | ENSG00000205231 |
        +--------------+--------------+-----------+-----------+--------------+-----------------+
        Stranded PyRanges object has 2,446 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.summary()
        +-------+------------------+--------------------+--------------------+-----------------------+
        |       |          pyrange |   coverage_forward |   coverage_reverse |   coverage_unstranded |
        |-------+------------------+--------------------+--------------------+-----------------------|
        | count |   2446           |               39   |               23   |                  32   |
        | mean  |   2291.92        |             7058.1 |            30078.6 |               27704.2 |
        | std   |  11906.9         |            10322.3 |            59467.7 |               67026.9 |
        | min   |      1           |               83   |              154   |                  83   |
        | 25%   |     90           |             1051   |             1204   |                1155   |
        | 50%   |    138           |             2541   |             6500   |                6343   |
        | 75%   |    382.25        |             7168   |            23778   |               20650.8 |
        | max   | 241726           |            43065   |           241726   |              291164   |
        | sum   |      5.60603e+06 |           275266   |           691807   |              886534   |
        +-------+------------------+--------------------+--------------------+-----------------------+

        >>> gr.summary(return_df=True, to_stdout=False)
                    pyrange  coverage_forward  coverage_reverse  coverage_unstranded
        count  2.446000e+03         39.000000         23.000000            32.000000
        mean   2.291918e+03       7058.102564      30078.565217         27704.187500
        std    1.190685e+04      10322.309347      59467.695265         67026.868647
        min    1.000000e+00         83.000000        154.000000            83.000000
        25%    9.000000e+01       1051.000000       1204.000000          1155.000000
        50%    1.380000e+02       2541.000000       6500.000000          6343.000000
        75%    3.822500e+02       7168.000000      23778.000000         20650.750000
        max    2.417260e+05      43065.000000     241726.000000        291164.000000
        sum    5.606031e+06     275266.000000     691807.000000        886534.000000
        """

        from pyranges.methods.summary import _summary

        return _summary(self, to_stdout, return_df)

    def tail(self, n: int = 8) -> "PyRanges":
        """Return the n last rows.

        Parameters
        ----------

        n : int, default 8

            Return n rows.

        Returns
        -------
        PyRanges

            PyRanges with the n last rows.

        See Also
        --------

        PyRanges.head : return the first rows
        PyRanges.sample : return random rows

        Examples
        --------

        >>> gr = pr.data.chipseq()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   | Start     | End       | Name       | Score     | Strand       |
        | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 212609534 | 212609559 | U0         | 0         | +            |
        | chr1         | 169887529 | 169887554 | U0         | 0         | +            |
        | chr1         | 216711011 | 216711036 | U0         | 0         | +            |
        | chr1         | 144227079 | 144227104 | U0         | 0         | +            |
        | ...          | ...       | ...       | ...        | ...       | ...          |
        | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
        | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
        | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
        | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.tail(3)
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chrY         |  13517892 |  13517917 | U0         |         0 | -            |
        | chrY         |   8010951 |   8010976 | U0         |         0 | -            |
        | chrY         |   7405376 |   7405401 | U0         |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        subsetter = np.zeros(len(self), dtype=np.bool_)
        subsetter[(len(self) - n) :] = True
        return self[subsetter]

    def tile(
        self, tile_size: int, overlap: bool = False, strand: Optional[bool] = None
    ) -> "PyRanges":
        """Return overlapping genomic tiles.

        The genome is divided into bookended tiles of length `tile_size` and one is returned per
        overlapping interval.

        Parameters
        ----------
        tile_size : int
            Length of the tiles.

        overlap : bool, default False

            Add column of nucleotide overlap to each tile.

        strand : bool, default None, i.e. auto

            Whether to do operations on chromosome/strand pairs or chromosomes. If None, will use
            chromosome/strand pairs if the PyRanges is stranded.

        **kwargs
            Additional keyword arguments to pass as keyword arguments to `f`

        Returns
        -------
        PyRanges

            Tiled PyRanges.

        See also
        --------

        pyranges.PyRanges.window : divide intervals into windows

        Examples
        --------

        >>> gr = pr.data.ensembl_gtf()[["Feature", "gene_name"]]
        >>> gr
        +--------------+--------------+-----------+-----------+--------------+-------------+
        | Chromosome   | Feature      | Start     | End       | Strand       | gene_name   |
        | (category)   | (category)   | (int64)   | (int64)   | (category)   | (object)    |
        |--------------+--------------+-----------+-----------+--------------+-------------|
        | 1            | gene         | 11868     | 14409     | +            | DDX11L1     |
        | 1            | transcript   | 11868     | 14409     | +            | DDX11L1     |
        | 1            | exon         | 11868     | 12227     | +            | DDX11L1     |
        | 1            | exon         | 12612     | 12721     | +            | DDX11L1     |
        | ...          | ...          | ...       | ...       | ...          | ...         |
        | 1            | gene         | 1173055   | 1179555   | -            | TTLL10-AS1  |
        | 1            | transcript   | 1173055   | 1179555   | -            | TTLL10-AS1  |
        | 1            | exon         | 1179364   | 1179555   | -            | TTLL10-AS1  |
        | 1            | exon         | 1173055   | 1176396   | -            | TTLL10-AS1  |
        +--------------+--------------+-----------+-----------+--------------+-------------+
        Stranded PyRanges object has 2,446 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.tile(200)
        +--------------+--------------+-----------+-----------+--------------+-------------+
        | Chromosome   | Feature      | Start     | End       | Strand       | gene_name   |
        | (category)   | (category)   | (int64)   | (int64)   | (category)   | (object)    |
        |--------------+--------------+-----------+-----------+--------------+-------------|
        | 1            | gene         | 11800     | 12000     | +            | DDX11L1     |
        | 1            | gene         | 12000     | 12200     | +            | DDX11L1     |
        | 1            | gene         | 12200     | 12400     | +            | DDX11L1     |
        | 1            | gene         | 12400     | 12600     | +            | DDX11L1     |
        | ...          | ...          | ...       | ...       | ...          | ...         |
        | 1            | exon         | 1175600   | 1175800   | -            | TTLL10-AS1  |
        | 1            | exon         | 1175800   | 1176000   | -            | TTLL10-AS1  |
        | 1            | exon         | 1176000   | 1176200   | -            | TTLL10-AS1  |
        | 1            | exon         | 1176200   | 1176400   | -            | TTLL10-AS1  |
        +--------------+--------------+-----------+-----------+--------------+-------------+
        Stranded PyRanges object has 30,538 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.tile(100, overlap=True)
        +--------------+--------------+-----------+-----------+--------------+-------------+---------------+
        | Chromosome   | Feature      | Start     | End       | Strand       | gene_name   | TileOverlap   |
        | (category)   | (category)   | (int64)   | (int64)   | (category)   | (object)    | (int64)       |
        |--------------+--------------+-----------+-----------+--------------+-------------+---------------|
        | 1            | gene         | 11800     | 11900     | +            | DDX11L1     | 32            |
        | 1            | gene         | 11900     | 12000     | +            | DDX11L1     | 100           |
        | 1            | gene         | 12000     | 12100     | +            | DDX11L1     | 100           |
        | 1            | gene         | 12100     | 12200     | +            | DDX11L1     | 100           |
        | ...          | ...          | ...       | ...       | ...          | ...         | ...           |
        | 1            | exon         | 1176000   | 1176100   | -            | TTLL10-AS1  | 100           |
        | 1            | exon         | 1176100   | 1176200   | -            | TTLL10-AS1  | 100           |
        | 1            | exon         | 1176200   | 1176300   | -            | TTLL10-AS1  | 100           |
        | 1            | exon         | 1176300   | 1176400   | -            | TTLL10-AS1  | 96            |
        +--------------+--------------+-----------+-----------+--------------+-------------+---------------+
        Stranded PyRanges object has 58,516 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        from pyranges.methods.windows import _tiles

        if strand is None:
            strand = self.strand_values_valid

        kwargs = {
            "strand": strand,
            "overlap": overlap,
            "sparse": {"self": False},
            "tile_size": tile_size,
        }

        return PyRanges(pyrange_apply_single(_tiles, self, **kwargs))

    def to_example(self, n: int = 10) -> Dict[str, List[Union[int, str]]]:
        """Return as dict.

        Used for easily creating examples for copy and pasting.

        Parameters
        ----------
        n : int, default 10
            Number of rows. Half is taken from the start, the other half from the end.

        See Also
        --------

        PyRanges.PyRanges : create PyRanges from dict

        Examples
        --------

        >>> gr = pr.data.chipseq()
        >>> gr
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   | Start     | End       | Name       | Score     | Strand       |
        | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 212609534 | 212609559 | U0         | 0         | +            |
        | chr1         | 169887529 | 169887554 | U0         | 0         | +            |
        | chr1         | 216711011 | 216711036 | U0         | 0         | +            |
        | chr1         | 144227079 | 144227104 | U0         | 0         | +            |
        | ...          | ...       | ...       | ...        | ...       | ...          |
        | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
        | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
        | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
        | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> d = gr.to_example(n=4)
        >>> d
        {'Chromosome': ['chr1', 'chr1', 'chrY', 'chrY'], 'Start': [212609534, 169887529, 8010951, 7405376], 'End': [212609559, 169887554, 8010976, 7405401], 'Name': ['U0', 'U0', 'U0', 'U0'], 'Score': [0, 0, 0, 0], 'Strand': ['+', '+', '-', '-']}
        >>> pr.PyRanges(d)
        +--------------+-----------+-----------+------------+-----------+--------------+
        | Chromosome   |     Start |       End | Name       |     Score | Strand       |
        | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
        |--------------+-----------+-----------+------------+-----------+--------------|
        | chr1         | 212609534 | 212609559 | U0         |         0 | +            |
        | chr1         | 169887529 | 169887554 | U0         |         0 | +            |
        | chrY         |   8010951 |   8010976 | U0         |         0 | -            |
        | chrY         |   7405376 |   7405401 | U0         |         0 | -            |
        +--------------+-----------+-----------+------------+-----------+--------------+
        Stranded PyRanges object has 4 rows and 6 columns from 2 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        nrows_half = int(min(n, len(self)) / 2)

        if n < len(self):
            first = self.head(nrows_half)
            last = self.tail(nrows_half)
            example = pr.concat([first, last])
        else:
            example = self

        d = {c: list(getattr(example, c)) for c in example.columns}

        return d

    def three_end(self) -> "PyRanges":
        """Return the 3'-end.

        The 3'-end is the start of intervals on the reverse strand and the end of intervals on the
        forward strand.

        Returns
        -------
        PyRanges
            PyRanges with the 3'.

        See Also
        --------
        PyRanges.five_end : return the five prime end

        Examples
        --------

        >>> d =  {'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...       'End': [5, 8], 'Strand': ['+', '-']}
        >>> gr = pr.PyRanges(d)
        >>> gr
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         5 | +            |
        | chr1         |         6 |         8 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.three_end()
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         4 |         5 | +            |
        | chr1         |         6 |         7 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        """

        assert self.strand_values_valid, "Need stranded pyrange to find 3'."
        kwargs = fill_kwargs({"strand": True})
        return pr.PyRanges(pyrange_apply_single(_tes, self, **kwargs))

    #     def to_bam(self, path=None, header=None, chromosome_sizes=None, chain=False):

    #         r"""Write to bam.

    #         Parameters
    #         ----------
    #         path : str, default None
    #             Where to write. If None, returns string representation.

    #         keep : bool, default True

    #             Whether to keep all columns, not just Chromosome, Start, End,
    #             Name, Score, Strand when writing.

    #         compression : str, compression type to use, by default infer based on extension.
    #             See pandas.DataFree.to_csv for more info.

    #         header : dict

    #             Header to use in the bamfile. See the pysam docs for how it should look.
    #             Or use the header attribute from another pyasam.AlignmentFile.

    #         chromosome_sizes : PyRanges or dict

    #             If dict: map of chromosome names to chromosome length.

    #         chain : bool, default False
    #             Whether to return the PyRanges after writing.

    #         Note
    #         ----

    #         The following pyranges columns are used when writing:

    #         Chromosome, Start, End, Strand, MapQ, Flag, QueryStart, QueryEnd, Name, Cigar, Quality

    #         Examples
    #         --------

    #         >>> header = {"SQ": [{"SN": 1, "LN": 249250621}]}

    #         >>> c = '''Name	Flag Chromosome Start End MapQ Cigar QuerySequence Quality
    # read1	115	1	142618765 142618790	255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:65536	ZL:i:25
    # read2	115	1	142618765 142618790  255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:214748	ZL:i:25
    # read3	115	1	142618765 142618790  255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:2147484	ZL:i:25
    # read4	115	1	142618765 142618790  255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:2147483647	ZL:i:25
    # read5	115	1	142618765 142618790  255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:-65536	ZL:i:25
    # read6	115	1	142618765 142618790  255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:-214748	ZL:i:25
    # read7	115	1	142618765 142618790  255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:-2147484	ZL:i:25
    # read8	115	1	142618765 142618790  255	25M	CGACCCACTCCGCCATTTTCATCCG	IIGIIIHIGIIFIIIIIIIGIGIII	NM:i:0ZP:i:-2147483647	ZL:i:25'''

    #         >>>
    #         """

    def to_bed(
        self, path: Optional[str] = None, keep: bool = True, compression: str = "infer"
    ) -> Optional[str]:
        r"""Write to bed.

        Parameters
        ----------
        path : str, default None
            Where to write. If None, returns string representation.

        keep : bool, default True

            Whether to keep all columns, not just Chromosome, Start, End,
            Name, Score, Strand when writing.

        compression : str, compression type to use, by default infer based on extension.
            See pandas.DataFree.to_csv for more info.

        Examples
        --------

        >>> d =  {'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...       'End': [5, 8], 'Strand': ['+', '-'], "Gene": [1, 2]}
        >>> gr = pr.PyRanges(d)
        >>> gr
        +--------------+-----------+-----------+--------------+-----------+
        | Chromosome   |     Start |       End | Strand       |      Gene |
        | (category)   |   (int64) |   (int64) | (category)   |   (int64) |
        |--------------+-----------+-----------+--------------+-----------|
        | chr1         |         1 |         5 | +            |         1 |
        | chr1         |         6 |         8 | -            |         2 |
        +--------------+-----------+-----------+--------------+-----------+
        Stranded PyRanges object has 2 rows and 5 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.to_bed()
        'chr1\t1\t5\t.\t.\t+\t1\nchr1\t6\t8\t.\t.\t-\t2\n'

        # File contents:
        chr1	1	5	.	.	+	1
        chr1	6	8	.	.	-	2

        Does not include noncanonical bed-column `Gene`:

        >>> gr.to_bed(keep=False)
        'chr1\t1\t5\t.\t.\t+\nchr1\t6\t8\t.\t.\t-\n'

        # File contents:
        chr1	1	5	.	.	+
        chr1	6	8	.	.	-

        >>> gr.to_bed("test.bed")

        >>> open("test.bed").readlines()
        ['chr1\t1\t5\t.\t.\t+\t1\n', 'chr1\t6\t8\t.\t.\t-\t2\n']
        """
        from pyranges.out import _to_bed

        result = _to_bed(self, path, keep=keep, compression=compression)

        return result

    def to_bigwig(
        self,
        path: None = None,
        chromosome_sizes: None = None,
        rpm: bool = True,
        divide: Optional[bool] = None,
        value_col: Optional[str] = None,
        dryrun: bool = False,
        chain: bool = False,
    ) -> Optional["PyRanges"]:
        """Write regular or value coverage to bigwig.

        Note
        ----

        To create one bigwig per strand, subset the PyRanges first.

        Parameters
        ----------
        path : str

            Where to write bigwig.

        chromosome_sizes : PyRanges or dict

            If dict: map of chromosome names to chromosome length.

        rpm : True

            Whether to normalize data by dividing by total number of intervals and multiplying by
            1e6.

        divide : bool, default False

            (Only useful with value_col) Divide value coverage by regular coverage and take log2.

        value_col : str, default None

            Name of column to compute coverage of.

        dryrun : bool, default False

            Return data that would be written without writing bigwigs.

        chain : bool, default False
            Whether to return the PyRanges after writing.

        Note
        ----

        Requires pybigwig to be installed.

        If you require more control over the normalization process, use pyranges.to_bigwig()

        See Also
        --------
        pyranges.to_bigwig : write pandas pd.DataFrame to bigwig.

        Examples
        --------

        >>> d =  {'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [1, 4, 6],
        ...       'End': [7, 8, 10], 'Strand': ['+', '-', '-'],
        ...       'Value': [10, 20, 30]}
        >>> gr = pr.PyRanges(d)
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

        >>> gr.to_bigwig(dryrun=True, rpm=False)
        +--------------+-----------+-----------+-------------+
        | Chromosome   |     Start |       End |       Score |
        | (category)   |   (int64) |   (int64) |   (float64) |
        |--------------+-----------+-----------+-------------|
        | chr1         |         1 |         4 |           1 |
        | chr1         |         4 |         6 |           2 |
        | chr1         |         6 |         7 |           3 |
        | chr1         |         7 |         8 |           2 |
        | chr1         |         8 |        10 |           1 |
        +--------------+-----------+-----------+-------------+
        Unstranded PyRanges object has 5 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.to_bigwig(dryrun=True, rpm=False, value_col="Value")
        +--------------+-----------+-----------+-------------+
        | Chromosome   |     Start |       End |       Score |
        | (category)   |   (int64) |   (int64) |   (float64) |
        |--------------+-----------+-----------+-------------|
        | chr1         |         1 |         4 |          10 |
        | chr1         |         4 |         6 |          30 |
        | chr1         |         6 |         7 |          60 |
        | chr1         |         7 |         8 |          50 |
        | chr1         |         8 |        10 |          30 |
        +--------------+-----------+-----------+-------------+
        Unstranded PyRanges object has 5 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.to_bigwig(dryrun=True, rpm=False, value_col="Value", divide=True)
        +--------------+-----------+-----------+-------------+
        | Chromosome   |     Start |       End |       Score |
        | (category)   |   (int64) |   (int64) |   (float64) |
        |--------------+-----------+-----------+-------------|
        | chr1         |         0 |         1 |   nan       |
        | chr1         |         1 |         4 |     3.32193 |
        | chr1         |         4 |         6 |     3.90689 |
        | chr1         |         6 |         7 |     4.32193 |
        | chr1         |         7 |         8 |     4.64386 |
        | chr1         |         8 |        10 |     4.90689 |
        +--------------+-----------+-----------+-------------+
        Unstranded PyRanges object has 6 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        from pyranges.out import _to_bigwig

        _chromosome_sizes = (
            pr.data.chromsizes() if chromosome_sizes is None else chromosome_sizes
        )

        result = _to_bigwig(
            self, path, _chromosome_sizes, rpm, divide, value_col, dryrun
        )

        if dryrun:
            return result

        if chain:
            return self
        else:
            return None

    def to_gff3(
        self,
        path: None = None,
        compression: str = "infer",
        map_cols: Optional[Dict[str, str]] = None,
    ) -> str | None:
        """Write to General Feature Format 3.

        The GFF format consists of a tab-separated file without header.
        GFF contains a fixed amount of columns, indicated below (names before ":").
        For each of these, PyRanges will use the corresponding column (names after ":").

        ``seqname: Chromosome
        source: Source
        type: Feature
        start: Start
        end: End
        score: Score
        strand: Strand
        phase: Frame
        attribute: autofilled``

        Columns which are not mapped to GFF columns are appended as a field
        in the attribute string (i.e. the last field).

        Parameters
        ----------
        path : str, default None, i.e. return string representation.

            Where to write file.

        compression : {‘infer’, ‘gzip’, ‘bz2’, ‘zip’, ‘xz’, None}, default "infer"

            Which compression to use. Uses file extension to infer by default.

        map_cols: dict, default None

            Override mapping between GFF and PyRanges fields for any number of columns.
            Format: ``{gff_column : pyranges_column}``
            If a mapping is found for the "attribute"` column, this not auto-filled

        Notes
        -----
        Nonexisting columns will be added with a '.' to represent the missing values.

        See Also
        --------
        pyranges.read_gff3 : read GFF3 files
        pyranges.to_gtf : write to GTF format

        Examples
        --------

        >>> d = {"Chromosome": [1] * 3, "Start": [1, 3, 5], "End": [4, 6, 9], "Feature": ["gene", "exon", "exon"]}
        >>> gr = pr.PyRanges(d)
        >>> gr.to_gff3()
        '1\\t.\\tgene\\t2\\t4\\t.\\t.\\t.\\t\\n1\\t.\\texon\\t4\\t6\\t.\\t.\\t.\\t\\n1\\t.\\texon\\t6\\t9\\t.\\t.\\t.\\t\\n'

        # How the file would look
        1	.	gene	2	4	.	.	.
        1	.	exon	4	6	.	.	.
        1	.	exon	6	9	.	.	.

        >>> gr.Gene = [1, 2, 3]
        >>> gr.function = ["a b", "c", "def"]
        >>> gr.to_gff3()
        '1\\t.\\tgene\\t2\\t4\\t.\\t.\\t.\\tGene=1;function=a b\\n1\\t.\\texon\\t4\\t6\\t.\\t.\\t.\\tGene=2;function=c\\n1\\t.\\texon\\t6\\t9\\t.\\t.\\t.\\tGene=3;function=def\\n'

        # How the file would look
        1	.	gene	2	4	.	.	.	Gene=1;function=a b
        1	.	exon	4	6	.	.	.	Gene=2;function=c
        1	.	exon	6	9	.	.	.	Gene=3;function=def

        >>> gr.the_frame = [0, 2, 1]
        >>> gr.tag = ['mRNA', 'CDS', 'CDS']
        >>> gr
        +--------------+-----------+-----------+------------+-----------+------------+-------------+------------+
        |   Chromosome |     Start |       End | Feature    |      Gene | function   |   the_frame | tag        |
        |   (category) |   (int64) |   (int64) | (object)   |   (int64) | (object)   |     (int64) | (object)   |
        |--------------+-----------+-----------+------------+-----------+------------+-------------+------------|
        |            1 |         1 |         4 | gene       |         1 | a b        |           0 | mRNA       |
        |            1 |         3 |         6 | exon       |         2 | c          |           2 | CDS        |
        |            1 |         5 |         9 | exon       |         3 | def        |           1 | CDS        |
        +--------------+-----------+-----------+------------+-----------+------------+-------------+------------+
        Unstranded PyRanges object has 3 rows and 8 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.

        >>> gr.to_gff3(map_cols={'phase':'the_frame', 'feature':'tag'})
        '1\\t.\\tmRNA\\t2\\t4\\t.\\t.\\t0\\tFeature=gene;Gene=1;function=a b\\n1\\t.\\tCDS\\t4\\t6\\t.\\t.\\t2\\tFeature=exon;Gene=2;function=c\\n1\\t.\\tCDS\\t6\\t9\\t.\\t.\\t1\\tFeature=exon;Gene=3;function=def\\n'

        # How the file would look
        1	.	mRNA	2	4	.	.	0	Gene=1;function=a b
        1	.	CDS	4	6	.	.	2	Gene=2;function=c
        1	.	CDS	6	9	.	.	1	Gene=3;function=def

        >>> gr.to_gff3(map_cols={'attribute':'Gene'})
        '1\\t.\\tgene\\t2\\t4\\t.\\t.\\t.\\tGene=1\\n1\\t.\\texon\\t4\\t6\\t.\\t.\\t.\\tGene=1\\n1\\t.\\texon\\t6\\t9\\t.\\t.\\t.\\tGene=1\\n'

        # How the file would look
        1	.	gene	2	4	.	.	.	Gene=1
        1	.	exon	4	6	.	.	.	Gene=1
        1	.	exon	6	9	.	.	.	Gene=1
        """

        from pyranges.out import _to_gff3

        return _to_gff3(self, path, compression=compression, map_cols=map_cols)

    def to_gtf(
        self,
        path: None = None,
        compression: str = "infer",
        map_cols: Optional[Dict[str, str]] = None,
    ) -> Optional[str]:
        """Write to Gene Transfer Format.

        The GTF format consists of a tab-separated file without header.
        It contains a fixed amount of columns, indicated below (names before ":").
        For each of these, PyRanges will use the corresponding column (names after ":").

        ``seqname: Chromosome
        source: Source
        type: Feature
        start: Start
        end: End
        score: Score
        strand: Strand
        frame: Frame
        attribute: auto-filled``

        Columns which are not mapped to GTF columns are appended as a field
        in the attribute string (i.e. the last field).

        Parameters
        ----------
        path : str, default None, i.e. return string representation.

            Where to write file.

        compression : {‘infer’, ‘gzip’, ‘bz2’, ‘zip’, ‘xz’, None}, default "infer"

            Which compression to use. Uses file extension to infer by default.

        map_cols: dict, default None

            Override mapping between GTF and PyRanges fields for any number of columns.
            Format: ``{gtf_column : pyranges_column}``
            If a mapping is found for the "attribute"` column, this not auto-filled

        Notes
        -----
        Nonexisting columns will be added with a '.' to represent the missing values.

        See Also
        --------
        pyranges.read_gtf : read GTF files
        pyranges.to_gff3 : write to GFF3 format

        Examples
        --------

        >>> d = {"Chromosome": [1] * 3, "Start": [1, 3, 5], "End": [4, 6, 9], "Feature": ["gene", "exon", "exon"]}
        >>> gr = pr.PyRanges(d)
        >>> gr.to_gtf()  # the raw string output
        '1\\t.\\tgene\\t2\\t4\\t.\\t.\\t.\\t\\n1\\t.\\texon\\t4\\t6\\t.\\t.\\t.\\t\\n1\\t.\\texon\\t6\\t9\\t.\\t.\\t.\\t\\n'

        # What the file contents look like:
        1	.	gene	2	4	.	.	.
        1	.	exon	4	6	.	.	.
        1	.	exon	6	9	.	.	.

        >>> gr.name = ["Tim", "Eric", "Endre"]
        >>> gr.prices = ["Cheap", "Premium", "Fine European"]
        >>> gr.to_gtf()  # the raw string output
        '1\\t.\\tgene\\t2\\t4\\t.\\t.\\t.\\tname "Tim"; prices "Cheap";\\n1\\t.\\texon\\t4\\t6\\t.\\t.\\t.\\tname "Eric"; prices "Premium";\\n1\\t.\\texon\\t6\\t9\\t.\\t.\\t.\\tname "Endre"; prices "Fine European";\\n'

        # What the file contents look like:
        1	.	gene	2	4	.	.	.	name "Tim"; prices "Cheap";
        1	.	exon	4	6	.	.	.	name "Eric"; prices "Premium";
        1	.	exon	6	9	.	.	.	name "Endre"; prices "Fine European";


        >>> gr.to_gtf(map_cols={"feature":"name", "attribute":"prices"})  # the raw string output
        '1\\t.\\tTim\\t2\\t4\\t.\\t.\\t.\\tprices "Cheap";\\n1\\t.\\tEric\\t4\\t6\\t.\\t.\\t.\\tprices "Premium";\\n1\\t.\\tEndre\\t6\\t9\\t.\\t.\\t.\\tprices "Fine European";\\n'

        # What the file contents look like:
        1	.	Tim	2	4	.	.	.	prices "Cheap";
        1	.	Eric	4	6	.	.	.	prices "Premium";
        1	.	Endre	6	9	.	.	.	prices "Fine European";
        """

        from pyranges.out import _to_gtf

        result = _to_gtf(self, path, compression=compression, map_cols=map_cols)

        return result

    def to_rle(
        self,
        value_col: Optional[str] = None,
        strand: Optional[bool] = None,
        rpm: bool = False,
    ) -> "RleDict":
        """Return as RleDict.

        Create collection of Rles representing the coverage or other numerical value.

        Parameters
        ----------
        value_col : str, default None
            Numerical column to create RleDict from.

        strand : bool, default None, i.e. auto
            Whether to treat strands serparately.

        rpm : bool, default False
            Normalize by multiplying with `1e6/(number_intervals)`.

        Returns
        -------
        pyrle.RleDict

            Rle with coverage or other info from the PyRanges.

        Examples
        --------

        >>> d = {'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 8, 5],
        ...      'End': [6, 9, 7], 'Score': [0.1, 5, 3.14], 'Strand': ['+', '+', '-']}
        >>> gr = pr.PyRanges(d)
        >>> gr.to_rle()
        chr1 +
        --
        +--------+-----+-----+-----+-----+
        | Runs   | 3   | 3   | 2   | 1   |
        |--------+-----+-----+-----+-----|
        | Values | 0.0 | 1.0 | 0.0 | 1.0 |
        +--------+-----+-----+-----+-----+
        Rle of length 9 containing 4 elements (avg. length 2.25)
        <BLANKLINE>
        chr1 -
        --
        +--------+-----+-----+
        | Runs   | 5   | 2   |
        |--------+-----+-----|
        | Values | 0.0 | 1.0 |
        +--------+-----+-----+
        Rle of length 7 containing 2 elements (avg. length 3.5)
        RleDict object with 2 chromosomes/strand pairs.

        >>> gr.to_rle(value_col="Score")
        chr1 +
        --
        +--------+-----+-----+-----+-----+
        | Runs   | 3   | 3   | 2   | 1   |
        |--------+-----+-----+-----+-----|
        | Values | 0.0 | 0.1 | 0.0 | 5.0 |
        +--------+-----+-----+-----+-----+
        Rle of length 9 containing 4 elements (avg. length 2.25)
        <BLANKLINE>
        chr1 -
        --
        +--------+-----+------+
        | Runs   | 5   | 2    |
        |--------+-----+------|
        | Values | 0.0 | 3.14 |
        +--------+-----+------+
        Rle of length 7 containing 2 elements (avg. length 3.5)
        RleDict object with 2 chromosomes/strand pairs.

        >>> gr.to_rle(value_col="Score", strand=False)
        chr1
        +--------+-----+-----+------+------+-----+-----+
        | Runs   | 3   | 2   | 1    | 1    | 1   | 1   |
        |--------+-----+-----+------+------+-----+-----|
        | Values | 0.0 | 0.1 | 3.24 | 3.14 | 0.0 | 5.0 |
        +--------+-----+-----+------+------+-----+-----+
        Rle of length 9 containing 6 elements (avg. length 1.5)
        Unstranded RleDict object with 1 chromosome.

        >>> gr.to_rle(rpm=True)
        chr1 +
        --
        +--------+-----+-------------------+-----+-------------------+
        | Runs   | 3   | 3                 | 2   | 1                 |
        |--------+-----+-------------------+-----+-------------------|
        | Values | 0.0 | 333333.3333333333 | 0.0 | 333333.3333333333 |
        +--------+-----+-------------------+-----+-------------------+
        Rle of length 9 containing 4 elements (avg. length 2.25)
        <BLANKLINE>
        chr1 -
        --
        +--------+-----+-------------------+
        | Runs   | 5   | 2                 |
        |--------+-----+-------------------|
        | Values | 0.0 | 333333.3333333333 |
        +--------+-----+-------------------+
        Rle of length 7 containing 2 elements (avg. length 3.5)
        RleDict object with 2 chromosomes/strand pairs.
        """

        if strand is None:
            strand = self.strand_values_valid

        from pyranges.methods.to_rle import _to_rle

        return _to_rle(self, value_col, strand=strand, rpm=rpm)

    def remove_strand(self) -> "PyRanges":
        """Remove strand.

        Note
        ----

        Removes Strand column even if PyRanges is not stranded.

        See Also
        --------

        PyRanges.stranded : whether PyRanges contains valid strand info.

        Examples
        --------

        >>> d =  {'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...       'End': [5, 8], 'Strand': ['+', '-']}
        >>> gr = pr.PyRanges(d)
        >>> gr
        +--------------+-----------+-----------+--------------+
        | Chromosome   |     Start |       End | Strand       |
        | (category)   |   (int64) |   (int64) | (category)   |
        |--------------+-----------+-----------+--------------|
        | chr1         |         1 |         5 | +            |
        | chr1         |         6 |         8 | -            |
        +--------------+-----------+-----------+--------------+
        Stranded PyRanges object has 2 rows and 4 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.remove_strand()
        +--------------+-----------+-----------+
        | Chromosome   |     Start |       End |
        | (category)   |   (int64) |   (int64) |
        |--------------+-----------+-----------|
        | chr1         |         1 |         5 |
        | chr1         |         6 |         8 |
        +--------------+-----------+-----------+
        Unstranded PyRanges object has 2 rows and 3 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome.
        """

        if "Strand" in self.columns:
            if (gr := self.drop(labels=[STRAND_COL], axis=1)) is not None:
                return PyRanges(gr)
        return self

    def window(self, window_size: int, strand: Optional[bool] = None) -> "PyRanges":
        """Return overlapping genomic windows.

        Windows of length `window_size` are returned.

        Parameters
        ----------
        window_size : int
            Length of the windows.

        strand : bool, default None, i.e. auto

            Whether to do operations on chromosome/strand pairs or chromosomes. If None, will use
            chromosome/strand pairs if the PyRanges is stranded.

        **kwargs
            Additional keyword arguments to pass as keyword arguments to `f`

        Returns
        -------
        PyRanges

            Tiled PyRanges.

        See also
        --------

        pyranges.PyRanges.tile : divide intervals into adjacent tiles.

        Examples
        --------

        >>> import pyranges as pr
        >>> gr = pr.PyRanges({"Chromosome": [1], "Start": [895], "End": [1259]})
        >>> gr
          Chromosome    Start      End
               int64    int64    int64
        ------------  -------  -------
                   1      895     1259


        >>> gr.window(200)
          Chromosome    Start      End
               int64    int64    int64
        ------------  -------  -------
                   1      895     1095
                   1     1095     1259


        >>> gr2 = pr.data.ensembl_gtf().get_with_loc_columns(["Feature", "gene_name"])
        >>> gr2
        Chromosome    Start    End      Strand      Feature     gene_name
        category      int64    int64    category    category    object
        ------------  -------  -------  ----------  ----------  -----------
        1             11868    14409    +           gene        DDX11L1
        1             11868    14409    +           transcript  DDX11L1
        1             11868    12227    +           exon        DDX11L1
        1             12612    12721    +           exon        DDX11L1
        ...           ...      ...      ...         ...         ...
        1             1173055  1179555  -           gene        TTLL10-AS1
        1             1173055  1179555  -           transcript  TTLL10-AS1
        1             1179364  1179555  -           exon        TTLL10-AS1
        1             1173055  1176396  -           exon        TTLL10-AS1




        >>> gr2 = pr.data.ensembl_gtf().get_with_loc_columns(["Feature", "gene_name"])
        >>> gr2.window(1000)
        +--------------+--------------+-----------+-----------+--------------+-------------+
        | Chromosome   | Feature      | Start     | End       | Strand       | gene_name   |
        | (category)   | (category)   | (int64)   | (int64)   | (category)   | (object)    |
        |--------------+--------------+-----------+-----------+--------------+-------------|
        | 1            | gene         | 11868     | 12868     | +            | DDX11L1     |
        | 1            | gene         | 12868     | 13868     | +            | DDX11L1     |
        | 1            | gene         | 13868     | 14409     | +            | DDX11L1     |
        | 1            | transcript   | 11868     | 12868     | +            | DDX11L1     |
        | ...          | ...          | ...       | ...       | ...          | ...         |
        | 1            | exon         | 1173055   | 1174055   | -            | TTLL10-AS1  |
        | 1            | exon         | 1174055   | 1175055   | -            | TTLL10-AS1  |
        | 1            | exon         | 1175055   | 1176055   | -            | TTLL10-AS1  |
        | 1            | exon         | 1176055   | 1176396   | -            | TTLL10-AS1  |
        +--------------+--------------+-----------+-----------+--------------+-------------+
        Stranded PyRanges object has 7,516 rows and 6 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        from pyranges.methods.windows import _windows

        if strand is None:
            strand = self.strand_values_valid

        kwargs = {
            "strand": strand,
            "sparse": {"self": False},
            "window_size": window_size,
        }

        df = pyrange_apply_single(_windows, self, **kwargs)
        return PyRanges(df)

    def get_with_loc_columns(self, key: str | Iterable[str]) -> "pr.PyRanges":
        """
        Return the requested columns and the genome location columns.

        Parameters
        ----------
        key
            Columns to return. Either a string or an iterable of strings.

        Returns
        -------
        PyRanges

            PyRanges with the requested columns.

        >>> gr = pr.PyRanges({"Chromosome": [1], "Start": [895], "End": [1259], "Strand": ["+"], "Score": [1], "Score2": [2]})
        >>> gr.get_with_loc_columns("Score2")
          Chromosome    Start      End  Strand      Score2
               int64    int64    int64  object       int64
        ------------  -------  -------  --------  --------
                   1      895     1259  +                2
        >>> gr.get_with_loc_columns(["Score2"])
          Chromosome    Start      End  Strand      Score2
               int64    int64    int64  object       int64
        ------------  -------  -------  --------  --------
                   1      895     1259  +                2
        """
        keys = [key] if isinstance(key, str) else key

        missing_loc_cols = [
            col
            for col in GENOME_LOC_COLS_WITH_STRAND
            if col not in keys and col in self.columns
        ]

        return PyRanges(super().__getitem__(missing_loc_cols + keys))

    def __getstate__(self):
        return self.dfs

    def __setstate__(self, d):
        self.__dict__["dfs"] = d

    @staticmethod
    def _zip_locationkey_and_data(
        keys: Iterable, dfs: Iterable[pd.DataFrame], strand: bool
    ) -> "PyRanges":
        """Zip keys and data into a PyRanges object.

        Helper method because MyPy has difficulty seeing that PyRanges keys are
        either list[str] or list[tuple[str, str]]. It considers them to be list[Union[str, tuple[str, str]]]
        which results in typecheck errors.
        """
        if strand:
            for k in keys:
                assert isinstance(k, tuple)
            return pr.PyRanges(dict(zip(keys, dfs)))
        else:
            for k in keys:
                assert isinstance(k, str)
            return pr.PyRanges(dict(zip(keys, dfs)))

    @property
    def _dfs_without_strand(self) -> Dict[str, pd.DataFrame]:
        """Return a dictionary of stranded dataframes."""
        assert not self.strand_values_valid, "PyRanges object is stranded"
        return {k: v for k, v in self.dfs.items() if isinstance(k, str)}

    @property
    def _dfs_with_strand(self) -> Dict[Tuple[str, str], pd.DataFrame]:
        """Return a dictionary of stranded dataframes."""
        assert self.strand_values_valid, "PyRanges object is not stranded"
        return {k: v for k, v in self.dfs.items() if isinstance(k, tuple)}


def _test():
    import doctest

    doctest.testmod()


if __name__ == "__main__":
    _test()
