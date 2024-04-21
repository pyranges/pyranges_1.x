"""Data structure for genomic intervals and their annotation."""

import logging
import sys
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional, cast

import numpy as np
import pandas as pd
from natsort import natsort, natsorted  # type: ignore[import]

import pyranges as pr
from pyranges.core.loci_getter import LociGetter
from pyranges.core.names import (
    CHROM_AND_STRAND_COLS,
    CHROM_COL,
    END_COL,
    FORWARD_STRAND,
    GENOME_LOC_COLS,
    GENOME_LOC_COLS_WITH_STRAND,
    JOIN_SUFFIX,
    NEAREST_DOWNSTREAM,
    NEAREST_UPSTREAM,
    OVERLAP_CONTAINMENT,
    OVERLAP_FIRST,
    PANDAS_COMPRESSION_TYPE,
    RANGE_COLS,
    SKIP_IF_EMPTY_LEFT,
    START_COL,
    STRAND_BEHAVIOR_AUTO,
    STRAND_BEHAVIOR_DEFAULT,
    STRAND_BEHAVIOR_OPPOSITE,
    STRAND_COL,
    TEMP_END_SLACK_COL,
    TEMP_NAME_COL,
    TEMP_NUM_COL,
    TEMP_START_SLACK_COL,
    TEMP_STRAND_COL,
    USE_STRAND_DEFAULT,
    VALID_BY_TYPES,
    VALID_COMBINE_OPTIONS,
    VALID_GENOMIC_STRAND_INFO,
    VALID_JOIN_TYPE,
    VALID_NEAREST_TYPE,
    VALID_OVERLAP_TYPE,
    VALID_STRAND_BEHAVIOR_TYPE,
    VALID_USE_STRAND_TYPE,
    BinaryOperation,
    CombineIntervalColumnsOperation,
    UnaryOperation,
)
from pyranges.core.parallelism import (
    _extend,
    _extend_grp,
    _tes,
    _tss,
)
from pyranges.core.pyranges_groupby import PyRangesDataFrameGroupBy
from pyranges.core.pyranges_helpers import (
    ensure_strand_behavior_options_valid,
    get_by_columns_including_chromosome_and_strand,
    group_keys_from_strand_behavior,
    group_keys_single,
    mypy_ensure_pyranges,
    strand_behavior_from_strand_and_validate,
    strand_from_strand_behavior,
    validate_and_convert_strand,
)
from pyranges.core.tostring import tostring
from pyranges.methods.merge import _merge
from pyranges.range_frame.range_frame import RangeFrame
from pyranges.range_frame.range_frame_validator import InvalidRangesReason

if TYPE_CHECKING:
    import pyfaidx  # type: ignore[import]
    from pyrle.rledict import Rledict  # type: ignore[import]


__all__ = ["PyRanges"]


logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class PyRanges(RangeFrame):
    """Two-dimensional representation of genomic intervals and their annotations.

    A PyRanges object must have the columns Chromosome, Start and End. These
    describe the genomic position and function as implicit row labels. A Strand
    column is optional and adds strand information to the intervals. Any other
    columns are allowed and are considered metadata.

    Operations between PyRanges align intervals based on their position.


    Parameters
    ----------
    You can initialize a PyRanges object like you would a pandas DataFrame, as long as the resulting DataFrame
    has the necessary columns (Chromosome, Start, End; Strand is optional).
    See examples below, and https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html for more information.

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
    index    |    Chromosome    Start      End
    int64    |    float64       float64    float64
    -------  ---  ------------  ---------  ---------
    PyRanges with 0 rows, 3 columns, and 1 index columns.
    Contains 0 chromosomes.

    You can initiatize PyRanges with a DataFrame:
    >>> df = pd.DataFrame({"Chromosome": ["chr1", "chr2"], "Start": [100, 200],
    ...                    "End": [150, 201]})
    >>> df
      Chromosome  Start  End
    0       chr1    100  150
    1       chr2    200  201
    >>> pr.PyRanges(df)
      index  |    Chromosome      Start      End
      int64  |    object          int64    int64
    -------  ---  ------------  -------  -------
          0  |    chr1              100      150
          1  |    chr2              200      201
    PyRanges with 2 rows, 3 columns, and 1 index columns.
    Contains 2 chromosomes.

    Or you can use a dictionary of iterables:
    >>> gr = pr.PyRanges({"Chromosome": [1, 1], "Strand": ["+", "-"], "Start": [1, 4], "End": [2, 27],
    ...                    "TP": [0, 1], "FP": [12, 11], "TN": [10, 9], "FN": [2, 3]})
    >>> gr
      index  |      Chromosome  Strand      Start      End       TP       FP       TN       FN
      int64  |           int64  object      int64    int64    int64    int64    int64    int64
    -------  ---  ------------  --------  -------  -------  -------  -------  -------  -------
          0  |               1  +               1        2        0       12       10        2
          1  |               1  -               4       27        1       11        9        3
    PyRanges with 2 rows, 8 columns, and 1 index columns.
    Contains 1 chromosomes and 2 strands.

    Operations that remove a column required for a PyRanges return a DataFrame instead
    >>> gr.drop("Chromosome", axis=1)
      Strand  Start  End  TP  FP  TN  FN
    0      +      1    2   0  12  10   2
    1      -      4   27   1  11   9   3

    >>> pr.PyRanges(dict(Chromosome=["chr1", "chr2"], Start=[1, 2], End=[2, 3]))
      index  |    Chromosome      Start      End
      int64  |    object          int64    int64
    -------  ---  ------------  -------  -------
          0  |    chr1                1        2
          1  |    chr2                2        3
    PyRanges with 2 rows, 3 columns, and 1 index columns.
    Contains 2 chromosomes.

    """

    def __new__(cls, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame":  # type: ignore[misc]
        """Create a new instance of a PyRanges object."""
        # __new__ is a special static method used for creating and
        # returning a new instance of a class. It is caladdled before
        # __init__ and is typically used in scenarios requiring
        # control over the creation of new instances

        # Logic to decide whether to return an instance of PyRanges or a DataFrame
        if not args and "data" not in kwargs:
            df = pd.DataFrame({k: [] for k in GENOME_LOC_COLS})
            df.__class__ = pr.PyRanges
            return df

        df = pd.DataFrame(kwargs.get("../data") or args[0])
        missing_any_required_columns = not set(GENOME_LOC_COLS).issubset({*df.columns})
        if missing_any_required_columns:
            return df

        return super().__new__(cls)

    def groupby(self, *args, **kwargs) -> "PyRangesDataFrameGroupBy":
        """Groupby PyRanges."""
        index_of_observed_in_args_list = 7
        if "observed" not in kwargs or len(args) < index_of_observed_in_args_list:
            kwargs["observed"] = True
        grouped = super().groupby(*args, **kwargs)
        return PyRangesDataFrameGroupBy(grouped)

    @property
    def _constructor(self) -> type:
        return pr.PyRanges

    @property
    def loci(self) -> LociGetter:
        """Get or set items in pyranges using .loci accessor.

        Notes
        -----
        In case of a 2-tuple with (val, slice), the method will first check if the val is in
        the chromosome col, and if so, it will subset on the matching rows. If val is not in
        the chromosome col it will look in the strand col.

        Parameters
        ----------
        key
            Genomic location (one or more of Chromosome, Strand, and Range),

        Returns
        -------
        PyRanges
            PyRanges with rows matching the location.

        Examples
        --------
        >>> import pyranges as pr
        >>> gr = pr.example_data.ensembl_gtf.get_with_loc_columns(["gene_id", "gene_name"])
        >>> gr
        index    |    Chromosome    Start    End      Strand      gene_id          gene_name
        int64    |    category      int64    int64    category    object           object
        -------  ---  ------------  -------  -------  ----------  ---------------  -----------
        0        |    1             11868    14409    +           ENSG00000223972  DDX11L1
        1        |    1             11868    14409    +           ENSG00000223972  DDX11L1
        2        |    1             11868    12227    +           ENSG00000223972  DDX11L1
        3        |    1             12612    12721    +           ENSG00000223972  DDX11L1
        ...      |    ...           ...      ...      ...         ...              ...
        7        |    1             120724   133723   -           ENSG00000238009  AL627309.1
        8        |    1             133373   133723   -           ENSG00000238009  AL627309.1
        9        |    1             129054   129223   -           ENSG00000238009  AL627309.1
        10       |    1             120873   120932   -           ENSG00000238009  AL627309.1
        PyRanges with 11 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.loci[1, "+", 12227:13000]
          index  |      Chromosome    Start      End  Strand      gene_id          gene_name
          int64  |        category    int64    int64  category    object           object
        -------  ---  ------------  -------  -------  ----------  ---------------  -----------
              0  |               1    11868    14409  +           ENSG00000223972  DDX11L1
              1  |               1    11868    14409  +           ENSG00000223972  DDX11L1
              3  |               1    12612    12721  +           ENSG00000223972  DDX11L1
        PyRanges with 3 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        >>> gr.loci[1, 14408:120000]
          index  |      Chromosome    Start      End  Strand      gene_id          gene_name
          int64  |        category    int64    int64  category    object           object
        -------  ---  ------------  -------  -------  ----------  ---------------  -----------
              0  |               1    11868    14409  +           ENSG00000223972  DDX11L1
              1  |               1    11868    14409  +           ENSG00000223972  DDX11L1
              4  |               1    13220    14409  +           ENSG00000223972  DDX11L1
              5  |               1   112699   112804  -           ENSG00000238009  AL627309.1
              6  |               1   110952   111357  -           ENSG00000238009  AL627309.1
        PyRanges with 5 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.loci[1, "-"]
          index  |      Chromosome    Start      End  Strand      gene_id          gene_name
          int64  |        category    int64    int64  category    object           object
        -------  ---  ------------  -------  -------  ----------  ---------------  -----------
              5  |               1   112699   112804  -           ENSG00000238009  AL627309.1
              6  |               1   110952   111357  -           ENSG00000238009  AL627309.1
              7  |               1   120724   133723  -           ENSG00000238009  AL627309.1
              8  |               1   133373   133723  -           ENSG00000238009  AL627309.1
              9  |               1   129054   129223  -           ENSG00000238009  AL627309.1
             10  |               1   120873   120932  -           ENSG00000238009  AL627309.1
        PyRanges with 6 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        >>> gr2 = pr.PyRanges({"Chromosome": ["chr1", "chr2"], "Start": [1, 2], "End": [4, 5], "Score": [9, 14], "Id": ["a", "b"]})
        >>> gr2.loci["chr2"]
          index  |    Chromosome      Start      End    Score  Id
          int64  |    object          int64    int64    int64  object
        -------  ---  ------------  -------  -------  -------  --------
              1  |    chr2                2        5       14  b
        PyRanges with 1 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr2.loci["chr2"] = gr2.loci["chr2"].assign(Chromosome="chr1")
        >>> gr2
          index  |    Chromosome      Start      End    Score  Id
          int64  |    object          int64    int64    int64  object
        -------  ---  ------------  -------  -------  -------  --------
              0  |    chr1                1        4        9  a
              1  |    chr1                2        5       14  b
        PyRanges with 2 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr2.loci["chr3"]
        Traceback (most recent call last):
        ...
        KeyError: 'Chromosome or strand "chr3" not found in PyRanges.'

        >>> gr2.loci["chr3", 1:2]
        Traceback (most recent call last):
        ...
        KeyError: 'Chromosome or strand chr3 not found in PyRanges.'
        >>> gr2.loci["Score"]
        Traceback (most recent call last):
        ...
        KeyError: 'Chromosome or strand "Score" not found in PyRanges.'

        >>> gr2.loci[["Score", "Id"]]
        Traceback (most recent call last):
        ...
        TypeError: The loci accessor does not accept a list. If you meant to retrieve columns, use gr.get_with_loc_columns instead.

        """
        return self._loci

    def sort_by_position(self) -> "pr.PyRanges":
        """Sort pyranges by Chromosome, Start, End, and possibly Strand.

        Examples
        --------
        >>> gr = pr.example_data.chipseq
        >>> gr
        index    |    Chromosome    Start      End        Name      Score    Strand
        int64    |    category      int64      int64      object    int64    category
        -------  ---  ------------  ---------  ---------  --------  -------  ----------
        0        |    chr8          28510032   28510057   U0        0        -
        1        |    chr7          107153363  107153388  U0        0        -
        2        |    chr5          135821802  135821827  U0        0        -
        3        |    chr14         19418999   19419024   U0        0        -
        ...      |    ...           ...        ...        ...       ...      ...
        16       |    chr9          120803448  120803473  U0        0        +
        17       |    chr6          89296757   89296782   U0        0        -
        18       |    chr1          194245558  194245583  U0        0        +
        19       |    chr8          57916061   57916086   U0        0        +
        PyRanges with 20 rows, 6 columns, and 1 index columns.
        Contains 15 chromosomes and 2 strands.

        >>> gr.sort_by_position()
        index    |    Chromosome    Start      End        Name      Score    Strand
        int64    |    category      int64      int64      object    int64    category
        -------  ---  ------------  ---------  ---------  --------  -------  ----------
        12       |    chr1          38457520   38457545   U0        0        +
        18       |    chr1          194245558  194245583  U0        0        +
        13       |    chr1          80668132   80668157   U0        0        -
        14       |    chr2          152562484  152562509  U0        0        -
        ...      |    ...           ...        ...        ...       ...      ...
        4        |    chr12         106679761  106679786  U0        0        -
        3        |    chr14         19418999   19419024   U0        0        -
        7        |    chr19         19571102   19571127   U0        0        +
        5        |    chr21         40099618   40099643   U0        0        +
        PyRanges with 20 rows, 6 columns, and 1 index columns.
        Contains 15 chromosomes and 2 strands.

        """
        self = mypy_ensure_pyranges(self.sort_values(([STRAND_COL] if self.has_strand else []) + RANGE_COLS))
        sorted_indexlike = np.array(
            natsort.order_by_index(
                [int(i) for i in self.index],
                natsort.index_natsorted(
                    self[CHROM_COL],
                ),
            ),
        )
        index = pd.Series(sorted_indexlike)
        return mypy_ensure_pyranges(
            self.reindex(
                index=index,
            ),
        )

    def __init__(self, *args, **kwargs) -> None:
        called_constructor_without_arguments = not args and "data" not in kwargs
        if called_constructor_without_arguments:
            # find out whether to include the strand column
            # also remove the strand key from kwargs since the constructor does not expect it.
            cols_to_use = GENOME_LOC_COLS_WITH_STRAND if kwargs.pop("strand", False) else GENOME_LOC_COLS
            # pass dict with the necessary pyranges cols to the constructor
            kwargs["data"] = {k: [] for k in cols_to_use}

        super().__init__(*args, **kwargs)

        self._loci = LociGetter(self)

    def _chrom_and_strand_info(self) -> str:
        num_chrom = self[CHROM_COL].nunique()
        strand_info = ""

        max_strands_to_show = 3

        if self.has_strand:
            num_strands = self[STRAND_COL].nunique()
            strands = f" and {num_strands} strands"
            if not self.strand_valid:
                nongenomic_strands = sorted(set(self[STRAND_COL]).difference(VALID_GENOMIC_STRAND_INFO))
                example_strands = ", ".join(
                    [str(s) for s in nongenomic_strands[:max_strands_to_show]]
                    + (["..."] if len(nongenomic_strands) > max_strands_to_show else []),
                )
                strands += f" (including non-genomic strands: {example_strands})"
            strand_info = strands

        return f"Contains {num_chrom} chromosomes{strand_info}"

    def __str__(self, **kwargs) -> str:
        r"""Return string representation.

        Examples
        --------
        >>> gr = pr.PyRanges({"Chromosome": [3, 2, 1], "Start": [0, 100, 250], "End": [10, 125, 251]})
        >>> gr
          index  |      Chromosome    Start      End
          int64  |           int64    int64    int64
        -------  ---  ------------  -------  -------
              0  |               3        0       10
              1  |               2      100      125
              2  |               1      250      251
        PyRanges with 3 rows, 3 columns, and 1 index columns.
        Contains 3 chromosomes.

        >>> gr.insert(3, "Strand", ["+", "-", "+"])
        >>> gr
          index  |      Chromosome    Start      End  Strand
          int64  |           int64    int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |               3        0       10  +
              1  |               2      100      125  -
              2  |               1      250      251  +
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        >>> gr.loc[:, "Strand"] = [".", "^", "/"]
        >>> gr
          index  |      Chromosome    Start      End  Strand
          int64  |           int64    int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |               3        0       10  .
              1  |               2      100      125  ^
              2  |               1      250      251  /
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 3 chromosomes and 3 strands (including non-genomic strands: ., /, ^).

        >>> gr2 = gr.copy()
        >>> gr2.loc[:, "Strand"] = ["+", "-", "X"]
        >>> gr = pr.concat([gr2, gr, gr])
        >>> gr  # If a PyRanges has more than eight rows the repr is truncated in the middle.
          index  |      Chromosome    Start      End  Strand
          int64  |           int64    int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |               3        0       10  +
              1  |               2      100      125  -
              2  |               1      250      251  X
              0  |               3        0       10  .
              1  |               2      100      125  ^
              2  |               1      250      251  /
              0  |               3        0       10  .
              1  |               2      100      125  ^
              2  |               1      250      251  /
        PyRanges with 9 rows, 4 columns, and 1 index columns (with 6 index duplicates).
        Contains 3 chromosomes and 6 strands (including non-genomic strands: ., /, X, ...).

        >>> gr = PyRanges({"Chromosome": [1], "Start": [1], "End": [2], "Strand": ["+"], "Name": ["Sonic The Hedgehog"], "gene_id": ["ENSG00000188976"], "short_gene_name": ["SHH"], "type": ["transcript"]})

        # The index is shown separated from the columns with |
        >>> gr.set_index(["gene_id"], append=True)
          level_0  gene_id          |      Chromosome    Start      End  Strand    Name                short_gene_name    ...
            int64  object           |           int64    int64    int64  object    object              object             ...
        ---------  ---------------  ---  ------------  -------  -------  --------  ------------------  -----------------  -----
                0  ENSG00000188976  |               1        1        2  +         Sonic The Hedgehog  SHH                ...
        PyRanges with 1 rows, 7 columns, and 2 index columns. (1 columns not shown: "type").
        Contains 1 chromosomes and 1 strands.

        >>> str_repr = gr.__str__(max_col_width=10, max_total_width=80)
        >>> print(str_repr)  # using print to show \\n as actual newlines
          index  |      Chromosome    Start      End  Strand    Name        gene_id     ...
          int64  |           int64    int64    int64  object    object      object      ...
        -------  ---  ------------  -------  -------  --------  ----------  ----------  -----
              0  |               1        1        2  +         Sonic T...  ENSG000...  ...
        PyRanges with 1 rows, 8 columns, and 1 index columns. (2 columns not shown: "short_gene_name", "type").
        Contains 1 chromosomes and 1 strands.

        >>> gr.insert(gr.shape[1], "Score", int(1e100))
        >>> gr.insert(gr.shape[1], "Score2", int(1e100))
        >>> str_repr = gr.__str__(max_col_width=10, max_total_width=80)
        >>> print(str_repr)
          index  |      Chromosome    Start      End  Strand    Name        gene_id     ...
          int64  |           int64    int64    int64  object    object      object      ...
        -------  ---  ------------  -------  -------  --------  ----------  ----------  -----
              0  |               1        1        2  +         Sonic T...  ENSG000...  ...
        PyRanges with 1 rows, 10 columns, and 1 index columns. (4 columns not shown: "short_gene_name", "type", "Score", ...).
        Contains 1 chromosomes and 1 strands.

        """
        str_repr = tostring(
            self,
            max_col_width=kwargs.get("max_col_width"),
            max_total_width=kwargs.get("max_total_width"),
        )
        str_repr = f"{str_repr}\n{self._chrom_and_strand_info()}."
        if reasons := InvalidRangesReason.formatted_reasons_list(self):
            str_repr = f"{str_repr}\nInvalid ranges:\n{reasons}"
        return str_repr

    def apply_single(
        self,
        function: UnaryOperation,
        by: VALID_BY_TYPES,
        use_strand: VALID_USE_STRAND_TYPE = USE_STRAND_DEFAULT,
        *,
        preserve_index: bool = False,
        **kwargs: Any,
    ) -> "pr.PyRanges":
        """Apply function to each group of intervals, defined by chromosome and optionally strand.

        Parameters
        ----------
        use_strand: {"auto", True, False}, default: "auto"
            Whether to use strand information when grouping.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        function : Callable
            Function that takes a PyRanges and optionally kwargs and returns a PyRanges.
            The function must accept a **kwargs argument. It may be used to extract useful information:
            use_strand = kwargs.get("use_strand", False)
            group = kwargs.get("__by__", {})
            # e.g. chromosome = group.get("Chromosome", None)
            # e.g. strand = group.get("Strand", "+")

        by : str or list of str or None
            Columns - in addition to chromosome and strand - to group by.

        preserve_index: bool
            Keep the old index. Only valid if the function preserves the index columns.

        kwargs : dict
            Arguments passed along to the function.

        """
        strand = validate_and_convert_strand(self, use_strand=use_strand)

        by = get_by_columns_including_chromosome_and_strand(self, by=by, use_strand=strand)
        return mypy_ensure_pyranges(
            super().apply_single(
                function=function,
                by=by,
                preserve_index=preserve_index,
                use_strand=use_strand,
                **kwargs,
            ),
        )

    def apply_pair(  # type: ignore[override]
        self: "PyRanges",
        other: "PyRanges",
        function: BinaryOperation,
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        by: VALID_BY_TYPES = None,
        **kwargs,
    ) -> "pr.PyRanges":
        """Apply function to pairs of overlapping intervals, by chromosome and optionally strand.

        Parameters
        ----------
        other : PyRanges
            Second PyRanges to apply function to.

        function : Callable
            Function that takes two PyRanges  and returns a PyRanges.
            The function shouldb accept a **kwargs argument. It may be used to extract useful information:
            group = kwargs.get("__by__", {})
            # e.g. chromosome = group.get("Chromosome", None)
            # e.g. strand = group.get("Strand", "+")

        strand_behavior : {"auto", "same", "opposite", "ignore"}, default "auto"
            Whether to consider overlaps of intervals on the same strand, the opposite or ignore strand
            information. The default, "auto", means use "same" if both PyRanges are stranded (see .strand_valid)
            otherwise ignore the strand information.

        by : str or list of str or None
            Additional columns - in addition to chromosome and strand - to group by.

        kwargs : dict
            Other arguments passed along to the function.

        """
        ensure_strand_behavior_options_valid(self, other, strand_behavior)
        by = self._by_to_list(by)

        if strand_behavior == STRAND_BEHAVIOR_OPPOSITE:
            self[TEMP_STRAND_COL] = self[STRAND_COL].replace({"+": "-", "-": "+"})
            other[TEMP_STRAND_COL] = other[STRAND_COL]

        grpby_ks = group_keys_from_strand_behavior(self, other, strand_behavior, by=by)

        res = mypy_ensure_pyranges(super().apply_pair(other, function, by=grpby_ks, **kwargs))

        if strand_behavior == STRAND_BEHAVIOR_OPPOSITE:
            res = res.drop_and_return(TEMP_STRAND_COL, axis="columns")

        return res

    def boundaries(
        self,
        transcript_id: str | list[str] | None = None,
        agg: dict[str, str | Callable] | None = None,
    ) -> "PyRanges":
        """Return the boundaries of groups of intervals (e.g. transcripts).

        While designed with transcript in mind, it can be used for any group of intervals.

        Parameters
        ----------
        transcript_id : str or list of str or None
            Name(s) of column(s) to group intervals.
            If None, intervals are grouped by chromosome, and strand if present and valid (see .strand_valid).

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
        >>> gr = pr.PyRanges({"Chromosome": [1, 1, 1], "Start": [1, 60, 110], "End": [40, 68, 130],
        ...                   "transcript_id": ["tr1", "tr1", "tr2"], "meta": ["a", "b", "c"]})
        >>>
        >>> gr["Length"] = gr.lengths()
        >>> gr
          index  |      Chromosome    Start      End  transcript_id    meta        Length
          int64  |           int64    int64    int64  object           object       int64
        -------  ---  ------------  -------  -------  ---------------  --------  --------
              0  |               1        1       40  tr1              a               39
              1  |               1       60       68  tr1              b                8
              2  |               1      110      130  tr2              c               20
        PyRanges with 3 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.boundaries("transcript_id")
          index  |      Chromosome    Start      End  transcript_id
          int64  |           int64    int64    int64  object
        -------  ---  ------------  -------  -------  ---------------
              0  |               1        1       68  tr1
              1  |               1      110      130  tr2
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.boundaries()
          index  |      Chromosome    Start      End
          int64  |           int64    int64    int64
        -------  ---  ------------  -------  -------
              0  |               1        1      130
        PyRanges with 1 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.boundaries("transcript_id", agg={"Length":"sum", "meta": ",".join})
          index  |      Chromosome    Start      End  transcript_id    meta        Length
          int64  |           int64    int64    int64  object           object       int64
        -------  ---  ------------  -------  -------  ---------------  --------  --------
              0  |               1        1       68  tr1              a,b             47
              1  |               1      110      130  tr2              c               20
        PyRanges with 2 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        from pyranges.methods.boundaries import _bounds

        # may be optimized: no need to split by chromosome/strands
        return self.apply_single(
            _bounds,
            by=self._by_to_list(transcript_id),
            agg=agg,
            use_strand=self.strand_valid,
        )

    @property
    def chromosomes(self) -> list[str]:
        """Return the list of unique chromosomes in this PyRanges, in natsorted order (e.g. chr2 < chr11)."""
        return natsorted(self[CHROM_COL].drop_duplicates())

    @property
    def chromosomes_and_strands(self) -> list[tuple[str, str]]:
        """Return the list of unique (chromosome, strand) pairs in this PyRanges in natsorted order (e.g. chr2 < chr11).

        Examples
        --------
        >>> gr = pr.PyRanges({"Chromosome": [1, 2, 2, 3], "Start": [1, 2, 3, 9], "End": [3, 3, 10, 12], "Strand": ["+", "-", "+", "-"]})
        >>> gr.chromosomes_and_strands
        [(1, '+'), (2, '+'), (2, '-'), (3, '-')]
        >>> gr.remove_strand().chromosomes_and_strands
        Traceback (most recent call last):
        ...
        ValueError: PyRanges has no strand column.

        """
        self._assert_strand_values_valid()
        return natsorted({*zip(self["Chromosome"], self["Strand"], strict=True)})

    def _assert_has_strand(self) -> None:
        if not self.has_strand:
            msg = "PyRanges has no strand column."
            raise ValueError(msg)

    def _assert_strand_values_valid(self) -> None:
        self._assert_has_strand()
        if not self.strand_valid:
            msg = f"PyRanges contains non-genomic strands. Only {VALID_GENOMIC_STRAND_INFO} are valid."
            raise ValueError(msg)

    def cluster(
        self,
        use_strand: VALID_USE_STRAND_TYPE = "auto",
        *,
        match_by: VALID_BY_TYPES = None,
        slack: int = 0,
        cluster_column: str = "Cluster",
        count_column: str | None = None,
    ) -> "PyRanges":
        """Give overlapping intervals a common id.

        Parameters
        ----------
        use_strand: {"auto", True, False}, default: "auto"
            Whether to cluster only intervals on the same strand.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        match_by : str or list, default None
            If provided, only intervals with an equal value in column(s) `match_by` may be considered as overlapping.

        slack : int, default 0
            Consider intervals separated by less than `slack` to be in the same cluster. If `slack`
            is negative, intervals overlapping less than `slack` are not considered to be in the
            same cluster.

        cluster_column:
            Name the cluster column added in output. Default: "Cluster"

        count_column:
            If a value is provided, add a column of counts with this name.

        Returns
        -------
        PyRanges
            PyRanges with an ID-column "Cluster" added.

        Warning
        -------

        Bookended intervals (i.e. the End of a PyRanges interval is the Start of
        another one) are by default considered to overlap.
        Avoid this with slack=-1.

        See Also
        --------
        PyRanges.merge: combine overlapping intervals into one

        Examples
        --------
        >>> gr = pr.PyRanges({"Chromosome": [1, 1, 2, 1, 1], "Start": [1, 2, 0, 3, 9],
        ...                    "End": [3, 3, 4, 10, 12], "Gene": [1, 2, 6, 3, 3]})
        >>> gr
          index  |      Chromosome    Start      End     Gene
          int64  |           int64    int64    int64    int64
        -------  ---  ------------  -------  -------  -------
              0  |               1        1        3        1
              1  |               1        2        3        2
              2  |               2        0        4        6
              3  |               1        3       10        3
              4  |               1        9       12        3
        PyRanges with 5 rows, 4 columns, and 1 index columns.
        Contains 2 chromosomes.

        >>> gr.cluster(cluster_column="ClusterId").sort_by_position()
          index  |      Chromosome    Start      End     Gene    ClusterId
          int64  |           int64    int64    int64    int64        int64
        -------  ---  ------------  -------  -------  -------  -----------
              0  |               1        1        3        1            0
              1  |               1        2        3        2            0
              3  |               1        3       10        3            0
              4  |               1        9       12        3            0
              2  |               2        0        4        6            1
        PyRanges with 5 rows, 5 columns, and 1 index columns.
        Contains 2 chromosomes.

        >>> gr.cluster(match_by=["Gene"], count_column="Counts")
          index  |      Chromosome    Start      End     Gene    Cluster    Counts
          int64  |           int64    int64    int64    int64      int64     int64
        -------  ---  ------------  -------  -------  -------  ---------  --------
              0  |               1        1        3        1          0         1
              1  |               1        2        3        2          1         1
              2  |               2        0        4        6          3         1
              3  |               1        3       10        3          2         2
              4  |               1        9       12        3          2         2
        PyRanges with 5 rows, 6 columns, and 1 index columns.
        Contains 2 chromosomes.

        """
        from pyranges.methods.cluster import _cluster

        strand = validate_and_convert_strand(self, use_strand=use_strand)
        _self = self.copy() if (not strand and self.has_strand) else self

        _by = [match_by] if isinstance(match_by, str) else ([*match_by] if match_by is not None else [])
        gr = _self.apply_single(
            _cluster,
            by=match_by,
            use_strand=strand,
            slack=slack,
            count_column=count_column,
            cluster_column=cluster_column,
            preserve_index=True,
        )
        gr[cluster_column] = gr.groupby(self.loc_columns + _by + [cluster_column]).ngroup()
        return mypy_ensure_pyranges(gr.reindex(self.index))

    def copy(self, *args, **kwargs) -> "pr.PyRanges":
        """Return a copy of the PyRanges."""
        return mypy_ensure_pyranges(super().copy(*args, **kwargs))

    def count_overlaps(
        self,
        other: "PyRanges",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        *,
        match_by: str | list[str] | None = None,
        overlap_col: str = "NumberOverlaps",
        keep_nonoverlapping: bool = True,
    ) -> "PyRanges":
        """Count number of overlaps per interval.

        For each interval in self, report how many intervals in 'other' overlap with it.

        Parameters
        ----------
        other: PyRanges
            Count overlaps with this PyRanges.

        match_by : str or list, default None
            If provided, only intervals with an equal value in column(s) `match_by` may be considered as overlapping.

        strand_behavior : {"auto", "same", "opposite", "ignore"}, default "auto"
            Whether to consider overlaps of intervals on the same strand, the opposite or ignore strand
            information. The default, "auto", means use "same" if both PyRanges are stranded (see .strand_valid)
            otherwise ignore the strand information.

        keep_nonoverlapping : bool, default True
            Keep intervals without overlaps.

        overlap_col : str, default "NumberOverlaps"
            Name of column with overlap counts.

        Returns
        -------
        PyRanges
            PyRanges with a column of overlaps added.

        See Also
        --------
        PyRanges.coverage: find coverage of PyRanges
        pyranges.count_overlaps: count overlaps from multiple PyRanges

        Examples
        --------
        >>> f1 = pr.example_data.f1.remove_nonloc_columns()
        >>> f1
          index  |    Chromosome      Start      End  Strand
          int64  |    category        int64    int64  category
        -------  ---  ------------  -------  -------  ----------
              0  |    chr1                3        6  +
              1  |    chr1                5        7  -
              2  |    chr1                8        9  +
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> f2 = pr.example_data.f2.remove_nonloc_columns()
        >>> f2
          index  |    Chromosome      Start      End  Strand
          int64  |    category        int64    int64  category
        -------  ---  ------------  -------  -------  ----------
              0  |    chr1                1        2  +
              1  |    chr1                6        7  -
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> f1.count_overlaps(f2, overlap_col="Count")
          index  |    Chromosome      Start      End  Strand        Count
          int64  |    category        int64    int64  category      int64
        -------  ---  ------------  -------  -------  ----------  -------
              0  |    chr1                3        6  +                 0
              2  |    chr1                8        9  +                 0
              1  |    chr1                5        7  -                 1
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        """
        from pyranges.methods.coverage import _number_overlapping

        return self.apply_pair(
            other,
            _number_overlapping,
            strand_behavior=strand_behavior,
            by=match_by,
            keep_nonoverlapping=keep_nonoverlapping,
            overlap_col=overlap_col,
            skip_if_empty=not keep_nonoverlapping,
        )

    def coverage(
        self,
        other: "PyRanges",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        *,
        overlap_col: str = "NumberOverlaps",
        fraction_col: str = "FractionOverlaps",
        match_by: VALID_BY_TYPES = None,
        keep_nonoverlapping: bool = True,
    ) -> "PyRanges":
        """Count number of overlaps and their fraction per interval.

        Count how many intervals in self overlap with those in other.

        Parameters
        ----------
        other: PyRanges
            Count overlaps from this PyRanges.

        match_by : str or list, default None
            If provided, only intervals with an equal value in column(s) `match_by` may be considered as overlapping.

        strand_behavior : {"auto", "same", "opposite", "ignore"}, default "auto"
            Whether to consider overlaps of intervals on the same strand, the opposite or ignore strand
            information. The default, "auto", means use "same" if both PyRanges are stranded (see .strand_valid)
            otherwise ignore the strand information.

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

        See Also
        --------
        pyranges.count_overlaps: count overlaps from multiple PyRanges

        Examples
        --------
        >>> f1 = pr.PyRanges({"Chromosome": [1, 1, 1], "Start": [3, 8, 5],
        ...                    "End": [6,  9, 7]})
        >>> f1
          index  |      Chromosome    Start      End
          int64  |           int64    int64    int64
        -------  ---  ------------  -------  -------
              0  |               1        3        6
              1  |               1        8        9
              2  |               1        5        7
        PyRanges with 3 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> f2 = pr.PyRanges({"Chromosome": [1, 1], "Start": [1, 6],
        ...                    "End": [2, 7]})
        >>> f2
          index  |      Chromosome    Start      End
          int64  |           int64    int64    int64
        -------  ---  ------------  -------  -------
              0  |               1        1        2
              1  |               1        6        7
        PyRanges with 2 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> f1.coverage(f2, overlap_col="C", fraction_col="F")
          index  |      Chromosome    Start      End        C          F
          int64  |           int64    int64    int64    int64    float64
        -------  ---  ------------  -------  -------  -------  ---------
              0  |               1        3        6        0        0
              1  |               1        8        9        0        0
              2  |               1        5        7        1        0.5
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        counts = self.count_overlaps(
            other,
            keep_nonoverlapping=True,
            overlap_col=overlap_col,
            strand_behavior=strand_behavior,
        )

        strand = strand_behavior != STRAND_BEHAVIOR_DEFAULT
        other = other.merge_overlaps(use_strand=strand, count_col="Count")

        from pyranges.methods.coverage import _coverage

        return counts.apply_pair(
            other,
            _coverage,
            strand_behavior=strand_behavior,
            by=match_by,
            fraction_col=fraction_col,
            keep_nonoverlapping=keep_nonoverlapping,
            overlap_col=overlap_col,
            skip_if_empty=not keep_nonoverlapping,
        )

    def extend(
        self,
        ext: int | None = None,
        ext_3: int | None = None,
        ext_5: int | None = None,
        transcript_id: str | list[str] | None = None,
    ) -> "PyRanges":
        """Extend the intervals from the ends.

        Parameters
        ----------
        ext: int or None
            Extend intervals by this amount from both ends.
        ext_3: int or None
            Extend intervals by this amount from the 3' end.
        ext_5: int or None
            Extend intervals by this amount from the 5' end.
        transcript_id : str or list of str, default: None
            group intervals into transcripts by these column name(s), so that the
            extension is applied only to the left-most and/or right-most interval.

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
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                3        6  +
              1  |    chr1                8        9  +
              2  |    chr1                5        7  -
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.


          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                3        7  +
              1  |    chr1                8       10  +
              2  |    chr1                4        7  -
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.extend(4)
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                0       10  +
              1  |    chr1                4       13  +
              2  |    chr1                1       11  -
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.


        >>> gr.extend(ext_3=1, ext_5=2)
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        7  +
              1  |    chr1                6       10  +
              2  |    chr1                4        9  -
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.extend(-1)
        Traceback (most recent call last):
        ...
        ValueError: Some intervals are negative or zero length after applying extend!

        >>> gr.extend(ext_3=1, ext_5=2)
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        7  +
              1  |    chr1                6       10  +
              2  |    chr1                4        9  -
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr['transcript_id']=['a', 'a', 'b']
        >>> gr.extend(transcript_id='transcript_id', ext_3=3)
          index  |    Chromosome      Start      End  Strand    transcript_id
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  ---------------
              0  |    chr1                3        6  +         a
              1  |    chr1                8       12  +         a
              2  |    chr1                2        7  -         b
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        """
        if (ext_3 or ext_5) and not self.strand_valid:
            msg = "PyRanges must be stranded to add 5/3-end specific extend."
            raise ValueError(msg)

        if ext is not None == (ext_3 is not None or ext_5 is not None):
            msg = "Must use at least one and not both of ext and ext3 or ext5."
            raise ValueError(msg)

        return (
            self.apply_single(
                _extend,
                by=group_keys_single(self, use_strand=self.strand_valid),
                ext=ext,
                ext_3=ext_3,
                ext_5=ext_5,
            )
            if not transcript_id
            else (
                self.apply_single(
                    _extend_grp,
                    by=group_keys_single(self, use_strand=self.strand_valid),
                    ext=ext,
                    ext_3=ext_3,
                    ext_5=ext_5,
                    group_by=transcript_id,
                )
            )
        )

    def five_end(
        self,
        transcript_id: str | list[str] | None = None,
    ) -> "PyRanges":
        """Return the five prime end of intervals.

        The five prime end is the start of a forward strand or the end of a reverse strand.

        Parameters
        ----------
        transcript_id : str or list of str, default: None
            Optional column name(s). If provided, the five prime end is calculated for each
            group of intervals.

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
        PyRanges.subsequence : return subintervals specified in relative genome-based coordinates
        PyRanges.spliced_subsequence : return subintervals specified in relative mRNA-based coordinates

        Examples
        --------
        >>> gr = pr.PyRanges({'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 10, 5], 'End': [9, 14, 7],
        ...                    'Strand': ["+", "+", "-"], 'Name': ['a', 'a', 'b']})
        >>> gr
          index  |    Chromosome      Start      End  Strand    Name
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  --------
              0  |    chr1                3        9  +         a
              1  |    chr1               10       14  +         a
              2  |    chr1                5        7  -         b
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.five_end()
          index  |    Chromosome      Start      End  Strand    Name
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  --------
              0  |    chr1                3        4  +         a
              1  |    chr1               10       11  +         a
              2  |    chr1                6        7  -         b
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.five_end(transcript_id='Name')
          index  |    Chromosome      Start      End  Strand    Name
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  --------
              0  |    chr1                3        4  +         a
              2  |    chr1                5        6  -         b
        PyRanges with 2 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        """
        if not self.strand_valid:
            msg = f"Need PyRanges with valid strands ({VALID_GENOMIC_STRAND_INFO}) to find 5'."
            raise AssertionError(msg)

        return (
            mypy_ensure_pyranges(self.apply_single(_tss, by=None, use_strand=True))
            if transcript_id is None
            else self.subsequence(transcript_id=transcript_id, start=0, end=1)
        )

    @property
    def loc_columns(self) -> list[str]:
        """Return the names of genomic location columns of this PyRanges (Chromosome, and Strand if present)."""
        return CHROM_AND_STRAND_COLS if self.has_strand else [CHROM_COL]

    @property
    def has_strand(self) -> bool:
        """Return whether PyRanges has a strand column.

        Does not check whether the strand column contains valid values.
        """
        return STRAND_COL in self.columns

    def join_ranges(
        self,
        other: "PyRanges",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        join_type: VALID_JOIN_TYPE = "inner",
        *,
        match_by: VALID_BY_TYPES = None,
        slack: int = 0,
        suffix: str = JOIN_SUFFIX,
        report_overlap: bool = False,
    ) -> "PyRanges":
        """Join PyRanges based on genomic overlap.

        Find pairs of overlapping intervals between two PyRanges (self and other) and combine their columns.
        Each row in the return PyRanges contains columns of both intervals, including their coordinates.
        By default, intervals without overlap are not reported.

        Parameters
        ----------
        other : PyRanges
            PyRanges to join.

        strand_behavior : {"auto", "same", "opposite", "ignore"}, default "auto"
            Whether to consider overlaps of intervals on the same strand, the opposite or ignore strand
            information. The default, "auto", means use "same" if both PyRanges are stranded (see .strand_valid)
            otherwise ignore the strand information.

        join_type : {"inner", "left", "right", "outer"}, default "inner"
            How to handle intervals without overlap. "inner" means only keep overlapping intervals.
            "left" keeps all intervals in self, "right" keeps all intervals in other, "outer" keeps both.
            For types other than "inner", intervals in self without overlaps will have NaN in columns from other,
            and/or vice versa.

        match_by : str or list, default None
            If provided, only intervals with an equal value in column(s) `match_by` may be joined.

        report_overlap : bool, default False
            Report amount of overlap in base pairs.

        slack : int, default 0
            Lengthen intervals in self before joining.

        suffix : str or tuple, default "_b"
            Suffix to give overlapping columns in other.

        apply_strand_suffix : bool, default None
            If first pyranges is unstranded, but the second is not, the first will be given a strand column.
            apply_strand_suffix makes the added strand column a regular data column instead by adding a suffix.

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

        See Also
        --------
        PyRanges.combine_interval_columns : give joined PyRanges new coordinates

        Examples
        --------
        >>> f1 = pr.PyRanges({'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 8, 5],
        ...                   'End': [6, 9, 7], 'Name': ['interval1', 'interval3', 'interval2']})
        >>> f1
          index  |    Chromosome      Start      End  Name
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  ---------
              0  |    chr1                3        6  interval1
              1  |    chr1                8        9  interval3
              2  |    chr1                5        7  interval2
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> f2 = pr.PyRanges({'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...                    'End': [2, 7], 'Name': ['a', 'b']})
        >>> f2
          index  |    Chromosome      Start      End  Name
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        2  a
              1  |    chr1                6        7  b
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> f1.join_ranges(f2)
          index  |    Chromosome      Start      End  Name       Chromosome_b      Start_b    End_b  Name_b
          int64  |    object          int64    int64  object     object              int64    int64  object
        -------  ---  ------------  -------  -------  ---------  --------------  ---------  -------  --------
              2  |    chr1                5        7  interval2  chr1                    6        7  b
        PyRanges with 1 rows, 8 columns, and 1 index columns.
        Contains 1 chromosomes.

        # Note that since some start and end columns are nan, a regular DataFrame is returned.

        >>> f1.join_ranges(f2, join_type="left")
          index  |    Chromosome      Start      End  Name       Chromosome_b      Start_b      End_b  Name_b
          int64  |    object          int64    int64  object     object            float64    float64  object
        -------  ---  ------------  -------  -------  ---------  --------------  ---------  ---------  --------
              2  |    chr1                5        7  interval2  chr1                    6          7  b
              0  |    chr1                3        6  interval1  nan                   nan        nan  nan
              1  |    chr1                8        9  interval3  nan                   nan        nan  nan
        PyRanges with 3 rows, 8 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> f1.join_ranges(f2, join_type="outer")
          index  |    Chromosome        Start        End  Name       Chromosome_b      Start_b      End_b  Name_b
          int64  |    object          float64    float64  object     object            float64    float64  object
        -------  ---  ------------  ---------  ---------  ---------  --------------  ---------  ---------  --------
              1  |    chr1                  5          7  interval2  chr1                    6          7  b
              0  |    chr1                  3          6  interval1  nan                   nan        nan  nan
              1  |    chr1                  8          9  interval3  nan                   nan        nan  nan
              0  |    nan                 nan        nan  nan        chr1                    1          2  a
        PyRanges with 4 rows, 8 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes.
        Invalid ranges:
          * 1 starts or ends are nan. See indexes: 0

        >>> gr = pr.PyRanges({"Chromosome": ["chr1", "chr2", "chr1", "chr3"], "Start": [1, 4, 10, 0],
        ...                   "End": [3, 9, 11, 1], "ID": ["a", "b", "c", "d"]})
        >>> gr
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        3  a
              1  |    chr2                4        9  b
              2  |    chr1               10       11  c
              3  |    chr3                0        1  d
        PyRanges with 4 rows, 4 columns, and 1 index columns.
        Contains 3 chromosomes.

        >>> gr2 = pr.PyRanges({"Chromosome": ["chr1", "chr1", "chr1"], "Start": [2, 2, 1], "End": [3, 9, 10],
        ...                     "ID": ["a", "b", "c"]})
        >>> gr2
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                2        3  a
              1  |    chr1                2        9  b
              2  |    chr1                1       10  c
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.join_ranges(gr2)
          index  |    Chromosome      Start      End  ID        Chromosome_b      Start_b    End_b  ID_b
          int64  |    object          int64    int64  object    object              int64    int64  object
        -------  ---  ------------  -------  -------  --------  --------------  ---------  -------  --------
              0  |    chr1                1        3  a         chr1                    1       10  c
              0  |    chr1                1        3  a         chr1                    2        9  b
              0  |    chr1                1        3  a         chr1                    2        3  a
        PyRanges with 3 rows, 8 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes.

        >>> gr.join_ranges(gr2, match_by="ID")
          index  |    Chromosome      Start      End  ID        Chromosome_b      Start_b    End_b  ID_b
          int64  |    object          int64    int64  object    object              int64    int64  object
        -------  ---  ------------  -------  -------  --------  --------------  ---------  -------  --------
              0  |    chr1                1        3  a         chr1                    2        3  a
        PyRanges with 1 rows, 8 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> bad = f1.join_ranges(f2, join_type="right")
        >>> bad
          index  |    Chromosome        Start        End  Name       Chromosome_b      Start_b    End_b  Name_b
          int64  |    object          float64    float64  object     object              int64    int64  object
        -------  ---  ------------  ---------  ---------  ---------  --------------  ---------  -------  --------
              1  |    chr1                  5          7  interval2  chr1                    6        7  b
              0  |    nan                 nan        nan  nan        chr1                    1        2  a
        PyRanges with 2 rows, 8 columns, and 1 index columns.
        Contains 1 chromosomes.
        Invalid ranges:
          * 1 starts or ends are nan. See indexes: 0
        >>> f2.join_ranges(bad)  # bad.join_ranges(f2) would not work either.
        Traceback (most recent call last):
        ...
        ValueError: Cannot perform function on invalid ranges (function was _both_dfs).

        With slack 1, bookended features are joined (see row 1):

        >>> f1.join_ranges(f2, slack=1)
          index  |    Chromosome      Start      End  Name       Chromosome_b      Start_b    End_b  Name_b
          int64  |    object          int64    int64  object     object              int64    int64  object
        -------  ---  ------------  -------  -------  ---------  --------------  ---------  -------  --------
              0  |    chr1                3        6  interval1  chr1                    6        7  b
              2  |    chr1                5        7  interval2  chr1                    6        7  b
        PyRanges with 2 rows, 8 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        from pyranges.methods.join import _both_dfs

        _self = self.copy()
        if slack:
            _self[TEMP_START_SLACK_COL] = _self.Start
            _self[TEMP_END_SLACK_COL] = _self.End

            _self = _self.extend(slack)

        gr: pd.DataFrame | PyRanges = _self.apply_pair(
            other,
            _both_dfs,
            strand_behavior=strand_behavior,
            by=match_by,
            join_type=join_type,
            report_overlap=report_overlap,
            suffix=suffix,
        )
        if slack and len(gr) > 0:
            gr[START_COL] = gr[TEMP_START_SLACK_COL]
            gr[END_COL] = gr[TEMP_END_SLACK_COL]
            gr = gr.drop_and_return([TEMP_START_SLACK_COL, TEMP_END_SLACK_COL], axis=1)

        return gr

    @property
    def length(self) -> int:
        """Return the total length of the intervals.

        See Also
        --------
        PyRanges.lengths : return the intervals lengths

        Examples
        --------
        >>> gr = pr.example_data.f1
        >>> gr
          index  |    Chromosome      Start      End  Name         Score  Strand
          int64  |    category        int64    int64  object       int64  category
        -------  ---  ------------  -------  -------  ---------  -------  ----------
              0  |    chr1                3        6  interval1        0  +
              1  |    chr1                5        7  interval2        0  -
              2  |    chr1                8        9  interval3        0  +
        PyRanges with 3 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.length
        6

        To find the length of the genome covered by the intervals, use merge first:

        >>> gr.merge_overlaps(use_strand=False).length
        5

        """
        lengths = self.lengths()
        length = lengths.sum()
        return int(length)

    def lengths(self) -> pd.Series:
        """Return the length of each interval.

        Returns
        -------
        pd.Series or dict of pd.Series with the lengths of each interval.

        See Also
        --------
        PyRanges.lengths : return the intervals lengths

        Examples
        --------
        >>> gr = pr.example_data.f1
        >>> gr
          index  |    Chromosome      Start      End  Name         Score  Strand
          int64  |    category        int64    int64  object       int64  category
        -------  ---  ------------  -------  -------  ---------  -------  ----------
              0  |    chr1                3        6  interval1        0  +
              1  |    chr1                5        7  interval2        0  -
              2  |    chr1                8        9  interval3        0  +
        PyRanges with 3 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.lengths()
        0    3
        1    2
        2    1
        dtype: int64

        >>> gr["Length"] = gr.lengths()
        >>> gr
          index  |    Chromosome      Start      End  Name         Score  Strand        Length
          int64  |    category        int64    int64  object       int64  category       int64
        -------  ---  ------------  -------  -------  ---------  -------  ----------  --------
              0  |    chr1                3        6  interval1        0  +                  3
              1  |    chr1                5        7  interval2        0  -                  2
              2  |    chr1                8        9  interval3        0  +                  1
        PyRanges with 3 rows, 7 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        """
        return self.End - self.Start

    def max_disjoint(
        self,
        use_strand: VALID_USE_STRAND_TYPE = "auto",
        *,
        slack: int = 0,
        match_by: VALID_BY_TYPES = None,
        **_,
    ) -> "PyRanges":
        """Find the maximal disjoint set of intervals.

        Returns a subset of the rows in self so that no two intervals overlap, choosing those that
        maximize the number of intervals in the result.

        Parameters
        ----------
        use_strand: {"auto", True, False}, default: "auto"
            Find the max disjoint set separately for each strand.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        slack : int, default 0
            Consider intervals within a distance of slack to be overlapping.

        match_by : str or list, default None
            If provided, only intervals with an equal value in column(s) `match_by` may be considered as overlapping.

        Returns
        -------
        PyRanges

            PyRanges with maximal disjoint set of intervals.

        See Also
        --------
        PyRanges.merge_overlaps : merge intervals into non-overlapping superintervals
        PyRanges.split : split intervals into non-overlapping subintervals

        Examples
        --------
        >>> gr = pr.example_data.f1
        >>> gr
          index  |    Chromosome      Start      End  Name         Score  Strand
          int64  |    category        int64    int64  object       int64  category
        -------  ---  ------------  -------  -------  ---------  -------  ----------
              0  |    chr1                3        6  interval1        0  +
              1  |    chr1                5        7  interval2        0  -
              2  |    chr1                8        9  interval3        0  +
        PyRanges with 3 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.max_disjoint(use_strand=False)
          index  |    Chromosome      Start      End  Name         Score  Strand
          int64  |    category        int64    int64  object       int64  category
        -------  ---  ------------  -------  -------  ---------  -------  ----------
              0  |    chr1                3        6  interval1        0  +
              2  |    chr1                8        9  interval3        0  +
        PyRanges with 2 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        """
        use_strand = validate_and_convert_strand(self, use_strand)
        from pyranges.methods.max_disjoint import _max_disjoint

        result = self.apply_single(_max_disjoint, by=match_by, use_strand=use_strand, preserve_index=True, slack=slack)

        # reordering as the original one
        common_index = self.index.intersection(result.index)
        result = result.reindex(common_index)

        return mypy_ensure_pyranges(result)

    def merge_overlaps(
        self,
        use_strand: VALID_USE_STRAND_TYPE = USE_STRAND_DEFAULT,
        *,
        count_col: str | None = None,
        match_by: VALID_BY_TYPES = None,
        slack: int = 0,
    ) -> "PyRanges":
        """Merge overlapping intervals into one.

        Parameters
        ----------
        use_strand: {"auto", True, False}, default: "auto"
            Only merge intervals on same strand.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        count : bool, default False
            Count intervals in each superinterval.

        count_col : str, default "Count"
            Name of column with counts.

        match_by : str or list, default None
            If provided, only intervals with an equal value in column(s) `match_by` may be considered as overlapping.

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
        >>> gr = pr.example_data.ensembl_gtf.get_with_loc_columns(["Feature", "gene_name"])
        >>> gr
        index    |    Chromosome    Start    End      Strand      Feature     gene_name
        int64    |    category      int64    int64    category    category    object
        -------  ---  ------------  -------  -------  ----------  ----------  -----------
        0        |    1             11868    14409    +           gene        DDX11L1
        1        |    1             11868    14409    +           transcript  DDX11L1
        2        |    1             11868    12227    +           exon        DDX11L1
        3        |    1             12612    12721    +           exon        DDX11L1
        ...      |    ...           ...      ...      ...         ...         ...
        7        |    1             120724   133723   -           transcript  AL627309.1
        8        |    1             133373   133723   -           exon        AL627309.1
        9        |    1             129054   129223   -           exon        AL627309.1
        10       |    1             120873   120932   -           exon        AL627309.1
        PyRanges with 11 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.merge_overlaps(count_col="Count")
          index  |      Chromosome    Start      End  Strand      Count
          int64  |          object    int64    int64  object      int64
        -------  ---  ------------  -------  -------  --------  -------
              0  |               1    11868    14409  +               5
              1  |               1   110952   111357  -               1
              2  |               1   112699   112804  -               1
              3  |               1   120724   133723  -               4
        PyRanges with 4 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.merge_overlaps(count_col="Count", match_by="gene_name")
          index  |      Chromosome    Start      End  Strand    gene_name      Count
          int64  |          object    int64    int64  object    object         int64
        -------  ---  ------------  -------  -------  --------  -----------  -------
              0  |               1    11868    14409  +         DDX11L1            5
              1  |               1   110952   111357  -         AL627309.1         1
              2  |               1   112699   112804  -         AL627309.1         1
              3  |               1   120724   133723  -         AL627309.1         4
        PyRanges with 4 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        """
        use_strand = validate_and_convert_strand(self, use_strand)
        _by = get_by_columns_including_chromosome_and_strand(self=self, by=match_by, use_strand=use_strand)
        return self.apply_single(_merge, by=_by, use_strand=use_strand, count_col=count_col, slack=slack)

    def nearest(
        self,
        other: "PyRanges",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "ignore",
        direction: VALID_NEAREST_TYPE = "any",
        *,
        suffix: str = JOIN_SUFFIX,
        overlap: bool = True,
    ) -> "PyRanges":
        """Find closest interval.

        For each interval in self PyRanges, the columns of the nearest interval in other PyRanges are appended.

        Parameters
        ----------
        other : PyRanges
            PyRanges to find nearest interval in.

        strand_behavior : {"auto", "same", "opposite", "ignore"}, default "auto"
            Whether to consider overlaps of intervals on the same strand, the opposite or ignore strand
            information. The default, "auto", means use "same" if both PyRanges are stranded (see .strand_valid)
            otherwise ignore the strand information.

        overlap : bool, default True
            Whether to include overlaps.

        direction : {"any", "upstream", "downstream"}, default "any", i.e. both directions
            Whether to only look for nearest in one direction. Always with respect to the PyRanges
            it is called on.

        suffix : str, default "_b"
            Suffix to give columns with shared name in other.

        Returns
        -------
        PyRanges

            A PyRanges with columns representing nearest interval horizontally appended.

        See Also
        --------
        PyRanges.join_ranges : Has a slack argument to find intervals within a distance.

        Examples
        --------
        >>> f1 = pr.example_data.f1.remove_nonloc_columns()
        >>> f1
          index  |    Chromosome      Start      End  Strand
          int64  |    category        int64    int64  category
        -------  ---  ------------  -------  -------  ----------
              0  |    chr1                3        6  +
              1  |    chr1                5        7  -
              2  |    chr1                8        9  +
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> f2 = pr.example_data.f2.remove_nonloc_columns()
        >>> f2
          index  |    Chromosome      Start      End  Strand
          int64  |    category        int64    int64  category
        -------  ---  ------------  -------  -------  ----------
              0  |    chr1                1        2  +
              1  |    chr1                6        7  -
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> f1.nearest(f2)
          index  |    Chromosome      Start      End  Strand        Start_b    End_b  Strand_b      Distance
          int64  |    category        int64    int64  category        int64    int64  category         int64
        -------  ---  ------------  -------  -------  ----------  ---------  -------  ----------  ----------
              1  |    chr1                5        7  -                   6        7  -                    0
              0  |    chr1                3        6  +                   6        7  -                    1
              1  |    chr1                8        9  +                   6        7  -                    2
        PyRanges with 3 rows, 8 columns, and 1 index columns (with 1 index duplicates).
        Contains 1 chromosomes and 2 strands.

        >>> f1.nearest(f2, direction="upstream")
          index  |    Chromosome      Start      End  Strand        Start_b    End_b  Strand_b      Distance
          int64  |    category        int64    int64  category        int64    int64  category         int64
        -------  ---  ------------  -------  -------  ----------  ---------  -------  ----------  ----------
              1  |    chr1                5        7  -                   6        7  -                    0
              0  |    chr1                3        6  +                   1        2  +                    2
              1  |    chr1                8        9  +                   6        7  -                    2
        PyRanges with 3 rows, 8 columns, and 1 index columns (with 1 index duplicates).
        Contains 1 chromosomes and 2 strands.

        """
        from pyranges.methods.nearest import _nearest

        if direction in {NEAREST_UPSTREAM, NEAREST_DOWNSTREAM} and not other.strand_valid:
            msg = "If doing upstream or downstream nearest, other pyranges must be stranded"
            raise AssertionError(msg)

        return self.apply_pair(
            other,
            _nearest,
            strand_behavior=strand_behavior,
            how=direction,
            overlap=overlap,
            suffix=suffix,
        )

    def overlap(  # type: ignore[override]
        self,
        other: "PyRanges",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        *,
        match_by: VALID_BY_TYPES = None,
        invert: bool = False,
        contained: bool = False,
        **_,
    ) -> "PyRanges":
        """Return overlapping intervals.

        Returns the intervals in self which overlap with those in other.

        Parameters
        ----------
        other : PyRanges
            PyRanges to find overlaps with.

        strand_behavior : {"auto", "same", "opposite", "ignore"}, default "auto"
            Whether to consider overlaps of intervals on the same strand, the opposite or ignore strand
            information. The default, "auto", means use "same" if both PyRanges are stranded (see .strand_valid)
            otherwise ignore the strand information.

        invert : bool, default False
            Whether to return the intervals without overlaps.

        contained : bool, default False
            Whether to report only intervals that are entirely contained in an interval of 'other'.

        match_by : str or list, default None
            If provided, only overlapping intervals with an equal value in column(s) `match_by` are reported.

        Returns
        -------
        PyRanges

            A PyRanges with overlapping intervals.

        See Also
        --------
        PyRanges.intersect : report overlapping subintervals
        PyRanges.set_intersect : set-intersect PyRanges (e.g. merge then intersect)

        Examples
        --------
        >>> gr = pr.PyRanges({"Chromosome": ["chr1", "chr1", "chr2", "chr1", "chr3"], "Start": [1, 1, 4, 10, 0],
        ...                    "End": [3, 3, 9, 11, 1], "ID": ["A", "a", "b", "c", "d"]})
        >>> gr2 = pr.PyRanges({"Chromosome": ["chr1", "chr1", "chr2"], "Start": [2, 2, 1], "End": [3, 9, 10]})
        >>> gr.overlap(gr2)
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        3  A
              1  |    chr1                1        3  a
              2  |    chr2                4        9  b
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 2 chromosomes.

        >>> gr.overlap(gr2, contained=True)
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              2  |    chr2                4        9  b
        PyRanges with 1 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.overlap(gr2, invert=True)
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              3  |    chr1               10       11  c
              4  |    chr3                0        1  d
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 2 chromosomes.

        >>> gr3 = pr.PyRanges({"Chromosome": 1, "Start": [2, 4], "End": [3, 5], "Strand": ["+", "-"]})
        >>> gr3
          index  |      Chromosome    Start      End  Strand
          int64  |           int64    int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |               1        2        3  +
              1  |               1        4        5  -
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr4 = pr.PyRanges({"Chromosome": 1, "Start": [0], "End": [10], "Strand": ["-"]})
        >>> gr3.overlap(gr4, strand_behavior="opposite")
          index  |      Chromosome    Start      End  Strand
          int64  |           int64    int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |               1        2        3  +
        PyRanges with 1 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        """
        from pyranges.methods.overlap import _overlap

        how = OVERLAP_CONTAINMENT if contained else OVERLAP_FIRST

        return self.apply_pair(
            other,
            _overlap,
            strand_behavior=strand_behavior,
            by=match_by,
            invert=invert,
            skip_if_empty=False if invert else SKIP_IF_EMPTY_LEFT,
            how=how,
        )

    def set_intersect(
        self,
        other: "PyRanges",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        multiple: VALID_OVERLAP_TYPE = "all",
    ) -> "PyRanges":
        """Return set-theoretical intersection.

        Like intersect, but both PyRanges are merged first.

        Parameters
        ----------
        other : PyRanges
            PyRanges to set-intersect.

        strand_behavior : {"auto", "same", "opposite", "ignore"}, default "auto"
            Whether to consider overlaps of intervals on the same strand, the opposite or ignore strand
            information. The default, "auto", means use "same" if both PyRanges are stranded (see .strand_valid)
            otherwise ignore the strand information.

        multiple : {"all", "first", "last"}, default "all"
            What to report when multiple merged intervals in 'other' overlap with the same merged interval in self.
            The default "all" reports all overlapping subintervals, which will have duplicate indices.
            "first" reports only, for each merged self interval, the overlapping 'other' subinterval with smallest Start
            "last" reports only the overlapping subinterval with the biggest End in 'other'

        Returns
        -------
        PyRanges

            A PyRanges with overlapping subintervals.

        See Also
        --------
        PyRanges.intersect : find overlapping subintervals
        PyRanges.overlap : report overlapping intervals

        Examples
        --------
        >>> r1 = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [5, 20, 40],"End": [10, 30, 50], "ID": ["a", "b", "c"]})
        >>> r1
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                5       10  a
              1  |    chr1               20       30  b
              2  |    chr1               40       50  c
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.


        >>> r2 = pr.PyRanges({"Chromosome": ["chr1"] * 4, "Start": [7, 18, 25, 28], "End": [9, 22, 33, 32]})
        >>> r2
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                7        9
              1  |    chr1               18       22
              2  |    chr1               25       33
              3  |    chr1               28       32
        PyRanges with 4 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.
        >>> r1.set_intersect(r2)
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                7        9
              1  |    chr1               20       22
              2  |    chr1               25       30
        PyRanges with 3 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> r1.set_intersect(r2, multiple='first')
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                7        9
              1  |    chr1               20       22
        PyRanges with 2 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        from pyranges.methods.overlap import _intersect

        use_strand = strand_from_strand_behavior(self, other, strand_behavior)
        self_clusters = self.merge_overlaps(use_strand=use_strand and self.has_strand)
        other_clusters = other.merge_overlaps(use_strand=use_strand and other.has_strand)
        return mypy_ensure_pyranges(
            self_clusters.apply_pair(
                other_clusters,
                _intersect,
                strand_behavior=strand_behavior,
                how=multiple,
            ).reset_index(drop=True),
        )

    def set_union(self, other: "PyRanges", strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto") -> "PyRanges":
        """Return set-theoretical union.

        Parameters
        ----------
        other : PyRanges
            PyRanges to do union with.

        strand_behavior : {"auto", "same", "opposite", "ignore"}, default "auto"
            Whether to consider overlaps of intervals on the same strand, the opposite or ignore strand
            information. The default, "auto", means use "same" if both PyRanges are stranded (see .strand_valid)
            otherwise ignore the strand information.

        Returns
        -------
        PyRanges

            A PyRanges with the union of intervals.

        See Also
        --------
        PyRanges.set_intersect : set-theoretical intersection
        PyRanges.overlap : report overlapping intervals

        Examples
        --------
        >>> gr = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [1, 4, 10],
        ...                    "End": [3, 9, 11], "ID": ["a", "b", "c"]})
        >>> gr
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        3  a
              1  |    chr1                4        9  b
              2  |    chr1               10       11  c
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr2 = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [2, 2, 9], "End": [3, 9, 10]})
        >>> gr2
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                2        3
              1  |    chr1                2        9
              2  |    chr1                9       10
        PyRanges with 3 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.set_union(gr2)
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                1       11
        PyRanges with 1 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        if self.empty and other.empty:
            return mypy_ensure_pyranges(self.copy())

        strand = strand_from_strand_behavior(self, other, strand_behavior)

        if not strand:
            self = self.remove_strand()
            other = other.remove_strand()

        gr = pr.concat([self, other])

        return gr.merge_overlaps(use_strand=strand)

    def sort_ranges(
        self,
        by: str | Iterable[str] | None = None,
        use_strand: VALID_USE_STRAND_TYPE = "auto",
        *,
        sort_descending: str | Iterable[str] | None = None,
        natsorting: bool = False,
        reverse: bool = False,
    ) -> "PyRanges":
        """Sort PyRanges according to Chromosome, Strand (if present), Start, and End; or by the specified columns.

        If PyRanges is stranded and use_strand is True, intervals on the negative strand are sorted in descending
        order, and End is considered before Start. This is to have a 5' to 3' order.
        For uses not covered by this function, use  DataFrame.sort_values().

        Parameters
        ----------
        by : str or list of str, default None
            If provided, sorting occurs by Chromosome, Strand (if present), *by, Start, and End.
            To prioritize columns differently (e.g. Strand before Chromosome), explicitly provide all columns
            in the desired order as part of the 'by' argument.
            You can prepend any column name with '-' to reverse the order of sorting for that column.

        use_strand: {"auto", True, False}, default: "auto"
            Whether negative strand intervals should be sorted in descending order, meaning 5' to 3'.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        sort_descending : str or list of str, default None
            A column name or list of column names to sort in descending order, instead of ascending.
            These may include column names in the 'by' argument, or those implicitly included (e.g. Chromosome).

        natsorting : bool, default False
            Whether to use natural sorting for Chromosome column, so that e.g. chr2 < chr11. Slows down sorting.

        reverse : bool, default False
            Whether to reverse the sort order.

        Returns
        -------
        PyRanges

            Sorted PyRanges. The index is preserved. Use .reset_index(drop=True) to reset the index.

        Examples
        --------
        >>> p = pr.PyRanges({"Chromosome": ["chr1", "chr1", "chr1", "chr1", "chr2", "chr11", "chr11", "chr1"],
        ...                  "Strand": ["+", "+", "-", "-", "+", "+", "+",  "+"],
        ...                  "Start": [40, 1, 10, 70, 300, 140, 160, 90],
        ...                  "End": [60, 11, 25, 80, 400, 152, 190, 100],
        ...                  "transcript_id":["t3", "t3", "t2", "t2", "t4", "t5", "t5", "t1"]})
        >>> p
          index  |    Chromosome    Strand      Start      End  transcript_id
          int64  |    object        object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              0  |    chr1          +              40       60  t3
              1  |    chr1          +               1       11  t3
              2  |    chr1          -              10       25  t2
              3  |    chr1          -              70       80  t2
              4  |    chr2          +             300      400  t4
              5  |    chr11         +             140      152  t5
              6  |    chr11         +             160      190  t5
              7  |    chr1          +              90      100  t1
        PyRanges with 8 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        >>> p.sort_ranges()
          index  |    Chromosome    Strand      Start      End  transcript_id
          int64  |    object        object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              1  |    chr1          +               1       11  t3
              0  |    chr1          +              40       60  t3
              7  |    chr1          +              90      100  t1
              3  |    chr1          -              70       80  t2
              2  |    chr1          -              10       25  t2
              5  |    chr11         +             140      152  t5
              6  |    chr11         +             160      190  t5
              4  |    chr2          +             300      400  t4
        PyRanges with 8 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        Do not sort negative strand intervals in descending order:
        >>> p.sort_ranges(use_strand=False)
          index  |    Chromosome    Strand      Start      End  transcript_id
          int64  |    object        object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              1  |    chr1          +               1       11  t3
              0  |    chr1          +              40       60  t3
              7  |    chr1          +              90      100  t1
              2  |    chr1          -              10       25  t2
              3  |    chr1          -              70       80  t2
              5  |    chr11         +             140      152  t5
              6  |    chr11         +             160      190  t5
              4  |    chr2          +             300      400  t4
        PyRanges with 8 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        Sort chromosomes in natural order:
        >>> p.sort_ranges(natsorting=True)
          index  |    Chromosome    Strand      Start      End  transcript_id
          int64  |    object        object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              1  |    chr1          +               1       11  t3
              0  |    chr1          +              40       60  t3
              7  |    chr1          +              90      100  t1
              3  |    chr1          -              70       80  t2
              2  |    chr1          -              10       25  t2
              4  |    chr2          +             300      400  t4
              5  |    chr11         +             140      152  t5
              6  |    chr11         +             160      190  t5
        PyRanges with 8 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        Sort by 'transcript_id' before than by columns Start and End (but after Chromosome and Strand):
        >>> p.sort_ranges(by='transcript_id')
          index  |    Chromosome    Strand      Start      End  transcript_id
          int64  |    object        object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              7  |    chr1          +              90      100  t1
              1  |    chr1          +               1       11  t3
              0  |    chr1          +              40       60  t3
              3  |    chr1          -              70       80  t2
              2  |    chr1          -              10       25  t2
              5  |    chr11         +             140      152  t5
              6  |    chr11         +             160      190  t5
              4  |    chr2          +             300      400  t4
        PyRanges with 8 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        Sort by 'transcript_id' before than by columns Strand, Start and End:
        >>> p.sort_ranges(by=['transcript_id', 'Strand'])
          index  |    Chromosome    Strand      Start      End  transcript_id
          int64  |    object        object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              7  |    chr1          +              90      100  t1
              3  |    chr1          -              70       80  t2
              2  |    chr1          -              10       25  t2
              1  |    chr1          +               1       11  t3
              0  |    chr1          +              40       60  t3
              5  |    chr11         +             140      152  t5
              6  |    chr11         +             160      190  t5
              4  |    chr2          +             300      400  t4
        PyRanges with 8 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        Same as before, but 'transcript_id' is sorted in descending order:
        >>> p.sort_ranges(by=['transcript_id', 'Strand'], sort_descending='transcript_id')
          index  |    Chromosome    Strand      Start      End  transcript_id
          int64  |    object        object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              1  |    chr1          +               1       11  t3
              0  |    chr1          +              40       60  t3
              3  |    chr1          -              70       80  t2
              2  |    chr1          -              10       25  t2
              7  |    chr1          +              90      100  t1
              5  |    chr11         +             140      152  t5
              6  |    chr11         +             160      190  t5
              4  |    chr2          +             300      400  t4
        PyRanges with 8 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        """
        by = [] if by is None else [by] if isinstance(by, str) else list(by)
        sort_descending = (
            []
            if sort_descending is None
            else [sort_descending]
            if isinstance(sort_descending, str)
            else list(sort_descending)
        )

        use_strand = validate_and_convert_strand(self, use_strand)

        cols_to_sort_for = (
            ([CHROM_COL] if CHROM_COL not in by else [])
            + ([STRAND_COL] if STRAND_COL not in by and STRAND_COL in self.columns else [])
            + by
            + ([START_COL] if START_COL not in by else [])
            + ([END_COL] if END_COL not in by else [])
        )

        missing_sort_descending_cols = [col for col in sort_descending if col not in cols_to_sort_for]
        if missing_sort_descending_cols:
            msg = "Sort_descending arguments must be among column names used for sorting! Not found: " + ", ".join(
                missing_sort_descending_cols,
            )
            raise ValueError(msg)

        ascending = [col not in sort_descending for col in cols_to_sort_for]
        if reverse:
            ascending = [not asc for asc in ascending]

        z = self.copy()
        if natsorting:
            natsort_fn = natsort.natsort_keygen()
            z = z.assign(**{TEMP_NAME_COL: natsort_fn(self[CHROM_COL])})
            cols_to_sort_for = [c if c != CHROM_COL else TEMP_NAME_COL for c in cols_to_sort_for]

        if not use_strand:
            z = z.sort_values(cols_to_sort_for, ascending=ascending)
        else:
            mask = z["Strand"] == "-"
            initial_starts = z.loc[mask, "Start"].copy()
            initial_ends = z.loc[mask, "End"].copy()

            # Swapping Start and End for negative strand intervals, and multiplying by -1 to sort in descending order
            z.loc[mask, "Start"], z.loc[mask, "End"] = (
                z.loc[mask, "End"].to_numpy() * -1,
                z.loc[mask, "Start"].to_numpy() * -1,
            )
            z = z.sort_values(cols_to_sort_for, ascending=ascending)

            # Swapping back
            z.loc[mask, "Start"], z.loc[mask, "End"] = initial_starts, initial_ends

        return mypy_ensure_pyranges(z.drop(TEMP_NAME_COL, axis=1) if natsorting else z)

    def sort_by_5_prime_ascending_and_3_prime_descending(
        self,
        *,
        reverse: bool = False,
    ) -> "PyRanges":
        """Sort by 5' end ascending and 3' end descending.

        Parameters
        ----------
        reverse : bool, default False
            Whether to reverse the sort order.

        Returns
        -------
        PyRanges

            Sorted PyRanges

        Examples
        --------
        >>> p = pr.PyRanges({"Chromosome": [1, 1, 1, 1, 2, 1, 1],
        ...                  "Strand": ["+", "+", "-", "-", "+", "+", "+"],
        ...                  "Start": [40, 1, 10, 70, 300, 140, 160],
        ...                  "End": [60, 11, 25, 80, 400, 152, 190],
        ...                  "transcript_id":["t3", "t3", "t2", "t2", "t4", "t1", "t1"]})
        >>> p
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              0  |               1  +              40       60  t3
              1  |               1  +               1       11  t3
              2  |               1  -              10       25  t2
              3  |               1  -              70       80  t2
              4  |               2  +             300      400  t4
              5  |               1  +             140      152  t1
              6  |               1  +             160      190  t1
        PyRanges with 7 rows, 5 columns, and 1 index columns.
        Contains 2 chromosomes and 2 strands.

        >>> p.sort_by_5_prime_ascending_and_3_prime_descending()
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              1  |               1  +               1       11  t3
              0  |               1  +              40       60  t3
              5  |               1  +             140      152  t1
              6  |               1  +             160      190  t1
              3  |               1  -              70       80  t2
              2  |               1  -              10       25  t2
              4  |               2  +             300      400  t4
        PyRanges with 7 rows, 5 columns, and 1 index columns.
        Contains 2 chromosomes and 2 strands.

        """

        def _sort_by_5_prime_ascending_and_3_prime_descending(df: "pr.PyRanges", **_) -> pd.DataFrame:
            positive_by = RANGE_COLS
            negative_by = RANGE_COLS[::-1]
            return pd.concat(
                [
                    df[df.Strand == "+"].sort_values(by=positive_by, ascending=not reverse),
                    df[df.Strand == "-"].sort_values(by=negative_by, ascending=reverse),
                ],
            )

        return self.apply_single(
            _sort_by_5_prime_ascending_and_3_prime_descending,
            by=None,
            preserve_index=True,
            use_strand=True,
        )

    def spliced_subsequence(
        self,
        start: int = 0,
        end: int | None = None,
        transcript_id: VALID_BY_TYPES = None,
        use_strand: VALID_USE_STRAND_TYPE = "auto",
        **_,
    ) -> "PyRanges":
        """Get subsequences of the intervals, using coordinates mapping to spliced transcripts (without introns).

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

        transcript_id : list of str, default None
            intervals are grouped by this/these ID column(s) beforehand, e.g. exons belonging to same transcripts


        use_strand: {"auto", True, False}, default: "auto"
            Whether strand is considered when interpreting the start and end arguments of this function.
            If True, counting is from the 5' end (the leftmost coordinate for + strand and the rightmost for - strand).
            If False, all intervals are processed like they reside on the + strand.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        Returns
        -------
        PyRanges
            Subregion of self, subsequenced as specified by arguments

        Note
        ----
        If the request goes out of bounds (e.g. requesting 100 nts for a 90nt region), only the existing portion is returned

        See Also
        --------
        PyRanges.subsequence : analogous to this method, but input coordinates refer to the unspliced transcript

        Examples
        --------
        >>> p  = pr.PyRanges({"Chromosome": [1, 1, 2, 2, 3],
        ...                   "Strand": ["+", "+", "-", "-", "+"],
        ...                   "Start": [1, 40, 10, 70, 140],
        ...                   "End": [11, 60, 25, 80, 152],
        ...                   "transcript_id":["t1", "t1", "t2", "t2", "t3"] })
        >>> p
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              0  |               1  +               1       11  t1
              1  |               1  +              40       60  t1
              2  |               2  -              10       25  t2
              3  |               2  -              70       80  t2
              4  |               3  +             140      152  t3
        PyRanges with 5 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        # Get the first 15 nucleotides of *each spliced transcript*, grouping exons by transcript_id:
        >>> p.spliced_subsequence(0, 15, transcript_id='transcript_id')
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              0  |               1  +               1       11  t1
              1  |               1  +              40       45  t1
              2  |               2  -              20       25  t2
              3  |               2  -              70       80  t2
              4  |               3  +             140      152  t3
        PyRanges with 5 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        # Get the last 20 nucleotides of each spliced transcript:
        >>> p.spliced_subsequence(-20, transcript_id='transcript_id')
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              1  |               1  +              40       60  t1
              2  |               2  -              10       25  t2
              3  |               2  -              70       75  t2
              4  |               3  +             140      152  t3
        PyRanges with 4 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        # Get region from 25 to 60 of each spliced transcript, or their existing subportion:
        >>> p.spliced_subsequence(25, 60, transcript_id='transcript_id')
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              1  |               1  +              55       60  t1
        PyRanges with 1 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        # Get region of each spliced transcript which excludes their first and last 3 nucleotides:
        >>> p.spliced_subsequence(3, -3, transcript_id='transcript_id')
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              0  |               1  +               4       11  t1
              1  |               1  +              40       57  t1
              2  |               2  -              13       25  t2
              3  |               2  -              70       77  t2
              4  |               3  +             143      149  t3
        PyRanges with 5 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        """
        if transcript_id is None:
            # in this case, the results of spliced_subsequence and subsequence are identical,
            # so we can optimize just one of them methods
            return self.subsequence(start=start, end=end, use_strand=use_strand)

        from pyranges.methods.spliced_subsequence import _spliced_subseq

        use_strand = validate_and_convert_strand(self, use_strand)

        sorted_p = self.sort_by_5_prime_ascending_and_3_prime_descending() if use_strand else self.sort_by_position()

        result = sorted_p.apply_single(
            _spliced_subseq,
            by=transcript_id,
            use_strand=use_strand,
            start=start,
            end=end,
            preserve_index=True,
        )

        # reordering as the original one
        common_index = self.index.intersection(result.index)
        result = result.reindex(common_index)

        return mypy_ensure_pyranges(result)

    def split(
        self,
        use_strand: VALID_USE_STRAND_TYPE = "auto",
        *,
        match_by: VALID_BY_TYPES = None,
        between: bool = False,
    ) -> "PyRanges":
        """Split into non-overlapping intervals.

        The output does not contain overlapping intervals, but intervals that are adjacent are not merged.
        No columns other than Chromosome, Start, End, and Strand (if present) are output.

        Parameters
        ----------
        use_strand: {"auto", True, False}, default: "auto"
            Whether to split only intervals on the same strand.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        between : bool, default False
            Output also intervals corresponding to the gaps between the intervals in self.

        match_by : str or list, default None
            If provided, only intervals with an equal value in column(s) `match_by` may be split.

        Returns
        -------
        PyRanges

            PyRanges with intervals split at overlap points.

        See Also
        --------
        pyranges.multioverlap : find overlaps with multiple PyRanges

        Examples
        --------
        >>> gr = pr.PyRanges({'Chromosome': ['chr1', 'chr1', 'chr1', 'chr1'], 'Start': [3, 5, 5, 11],
        ...                   'End': [6, 9, 7, 12], 'Strand': ['+', '+', '-', '-']})
        >>> gr
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                3        6  +
              1  |    chr1                5        9  +
              2  |    chr1                5        7  -
              3  |    chr1               11       12  -
        PyRanges with 4 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.split()
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                3        5  +
              1  |    chr1                5        6  +
              2  |    chr1                6        9  +
              3  |    chr1                5        7  -
              5  |    chr1               11       12  -
        PyRanges with 5 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.split(between=True)
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                3        5  +
              1  |    chr1                5        6  +
              2  |    chr1                6        9  +
              3  |    chr1                5        7  -
              4  |    chr1                7       11  -
              5  |    chr1               11       12  -
        PyRanges with 6 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.split(use_strand=False)
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                3        5
              1  |    chr1                5        6
              2  |    chr1                6        7
              3  |    chr1                7        9
              5  |    chr1               11       12
        PyRanges with 5 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.split(use_strand=False, between=True)
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                3        5
              1  |    chr1                5        6
              2  |    chr1                6        7
              3  |    chr1                7        9
              4  |    chr1                9       11
              5  |    chr1               11       12
        PyRanges with 6 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr['ID'] = ['a', 'b', 'a', 'c']
        >>> gr
          index  |    Chromosome      Start      End  Strand    ID
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  --------
              0  |    chr1                3        6  +         a
              1  |    chr1                5        9  +         b
              2  |    chr1                5        7  -         a
              3  |    chr1               11       12  -         c
        PyRanges with 4 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.split(use_strand=False, match_by='ID')
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                3        5
              1  |    chr1                5        6
              2  |    chr1                6        7
              3  |    chr1                5        9
              4  |    chr1               11       12
        PyRanges with 5 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        from pyranges.methods.split import _split

        use_strand = validate_and_convert_strand(self, use_strand=use_strand)
        df = self.apply_single(
            _split,
            by=match_by,
            preserve_index=False,
            use_strand=use_strand,
        )
        if not between:
            df = df.overlap(
                self,
                how=OVERLAP_FIRST,
                strand_behavior=strand_behavior_from_strand_and_validate(self, use_strand),
            )

        return df

    @property
    def strand_valid(self) -> bool:
        """Whether PyRanges has valid strand info.

        Values other than '+' and '-' in the Strand column are not considered valid.
        A PyRanges without a Strand column is also not considered to have valid strand info.

        See Also
        --------
        PyRanges.has_strand : whether a Strand column is present
        PyRanges.make_strand_valid : make the strand information in PyRanges valid

        Examples
        --------
        >>> gr = pr.PyRanges({'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...                   'End': [5, 8], 'Strand': ['+', '.']})
        >>> gr
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        5  +
              1  |    chr1                6        8  .
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands (including non-genomic strands: .).
        >>> gr.strand_valid  # invalid strand value: '.'
        False

        >>> "Strand" in gr.columns
        True

        """
        if STRAND_COL not in self.columns and len(self) > 0:
            return False
        return bool(self[STRAND_COL].isin(VALID_GENOMIC_STRAND_INFO).all())

    def make_strand_valid(self) -> "PyRanges":
        """Make the strand information in PyRanges valid.

        Convert all invalid Strand values (those other than "+" and "-") to positive stranded values "+".
        If the Strand column is not present, add it with all values set to "+".

        Returns
        -------
        PyRanges
            PyRanges with valid strand information.

        See Also
        --------
        PyRanges.strand_valid : whether PyRanges has valid strand info
        PyRanges.remove_strand : remove the Strand column from PyRanges

        Examples
        --------
        >>> gr = pr.PyRanges({'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...                   'End': [5, 8], 'Strand': ['-', '.']})
        >>> gr
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        5  -
              1  |    chr1                6        8  .
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands (including non-genomic strands: .).

        >>> gr.strand_valid
        False

        >>> gr.make_strand_valid()
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        5  -
              1  |    chr1                6        8  +
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr2 = pr.PyRanges({'Chromosome': ['chr1', 'chr1'], 'Start': [5, 22],
        ...                    'End': [15, 30]})
        >>> gr2
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                5       15
              1  |    chr1               22       30
        PyRanges with 2 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr2.make_strand_valid()
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                5       15  +
              1  |    chr1               22       30  +
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        """
        if STRAND_COL not in self.columns:
            result = self.copy().assign(**{STRAND_COL: "+"})
        else:
            result = self.copy()
            result.loc[~(self[STRAND_COL].isin(VALID_GENOMIC_STRAND_INFO)), STRAND_COL] = "+"

        return mypy_ensure_pyranges(result)

    def subsequence(
        self,
        start: int = 0,
        end: int | None = None,
        transcript_id: VALID_BY_TYPES = None,
        use_strand: VALID_USE_STRAND_TYPE = "auto",
        **_,
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

        transcript_id : list of str, default None
            intervals are grouped by this/these ID column(s) beforehand, e.g. exons belonging to same transcripts

        use_strand: {"auto", True, False}, default: "auto"
            Whether strand is considered when interpreting the start and end arguments of this function.
            If True, counting is from the 5' end (the leftmost coordinate for + strand and the rightmost for - strand).
            If False, all intervals are processed like they reside on the + strand.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        Returns
        -------
        PyRanges
            Subregion of self, subsequenced as specified by arguments

        Note
        ----
        If the request goes out of bounds (e.g. requesting 100 nts for a 90nt region), only the existing portion is returned

        See Also
        --------
        PyRanges.spliced_subsequence : analogous to this method, but intronic regions are not counted, so that input coordinates refer to the spliced transcript


        Examples
        --------
        >>> p  = pr.PyRanges({"Chromosome": [1, 1, 2, 2, 3],
        ...                   "Strand": ["+", "+", "-", "-", "+"],
        ...                   "Start": [1, 40, 2, 30, 140],
        ...                   "End": [20, 60, 13, 45, 155],
        ...                   "transcript_id":["t1", "t1", "t2", "t2", "t3"] })
        >>> p
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              0  |               1  +               1       20  t1
              1  |               1  +              40       60  t1
              2  |               2  -               2       13  t2
              3  |               2  -              30       45  t2
              4  |               3  +             140      155  t3
        PyRanges with 5 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        # Get the first 10 nucleotides (at the 5') of *each interval* (each line of the dataframe):
        >>> p.subsequence(0, 10)
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              0  |               1  +               1       11  t1
              1  |               1  +              40       50  t1
              2  |               2  -               2       12  t2
              3  |               2  -              30       40  t2
              4  |               3  +             140      150  t3
        PyRanges with 5 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        # Get the first 10 nucleotides of *each transcript*, grouping exons by transcript_id:
        >>> p.subsequence(0, 10, transcript_id='transcript_id')
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              0  |               1  +               1       11  t1
              2  |               2  -               2       12  t2
              4  |               3  +             140      150  t3
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        # Get the last 20 nucleotides of each transcript:
        >>> p.subsequence(-20, transcript_id='transcript_id')
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              1  |               1  +              40       60  t1
              3  |               2  -              30       45  t2
              4  |               3  +             140      155  t3
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        # Get region from 30 to 330 of each transcript, or their existing subportion:
        >>> p.subsequence(30, 300, transcript_id='transcript_id')
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              1  |               1  +              40       60  t1
              3  |               2  -              32       45  t2
        PyRanges with 2 rows, 5 columns, and 1 index columns.
        Contains 2 chromosomes and 2 strands.

        """
        from pyranges.methods.subsequence import _subseq

        result = self.apply_single(
            _subseq,
            by=transcript_id,
            use_strand=use_strand,
            start=start,
            end=end,
            preserve_index=True,
        )

        return mypy_ensure_pyranges(result)

    def subtract_ranges(
        self,
        other: "pr.PyRanges",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        *,
        match_by: VALID_BY_TYPES = None,
    ) -> "pr.PyRanges":
        """Subtract intervals, i.e. return non-overlapping subintervals.

        Identify intervals in other that overlap with intervals in self; return self with the overlapping parts removed.

        Parameters
        ----------
        other:
            PyRanges to subtract.

        strand_behavior: "auto", "same", "opposite", "ignore"
            How to handle strand information. "auto" means use "same" if both PyRanges are stranded,
            otherwise ignore the strand information. "same" means only subtract intervals on the same strand.
            "opposite" means only subtract intervals on the opposite strand. "ignore" means ignore strand

        match_by : str or list, default None
            If provided, only intervals with an equal value in column(s) `match_by` may be considered as overlapping.

        See Also
        --------
        PyRanges.overlap : use with invert=True to return all intervals without overlap

        Examples
        --------
        >>> gr = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [1, 4, 10],
        ...                    "End": [3, 9, 11], "ID": ["a", "b", "c"]})
        >>> gr2 = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [2, 2, 9], "End": [3, 9, 10]})
        >>> gr
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        3  a
              1  |    chr1                4        9  b
              2  |    chr1               10       11  c
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr2
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                2        3
              1  |    chr1                2        9
              2  |    chr1                9       10
        PyRanges with 3 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.subtract_ranges(gr2)
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        2  a
              2  |    chr1               10       11  c
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr['tag'] = ['x', 'y', 'z']
        >>> gr2['tag'] = ['x', 'w', 'z']
        >>> gr
          index  |    Chromosome      Start      End  ID        tag
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  --------
              0  |    chr1                1        3  a         x
              1  |    chr1                4        9  b         y
              2  |    chr1               10       11  c         z
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr2
          index  |    Chromosome      Start      End  tag
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                2        3  x
              1  |    chr1                2        9  w
              2  |    chr1                9       10  z
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.subtract_ranges(gr2, match_by="tag")
          index  |    Chromosome      Start      End  ID        tag
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  --------
              0  |    chr1                1        2  a         x
              1  |    chr1                4        9  b         y
              2  |    chr1               10       11  c         z
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        from pyranges.methods.subtraction import _subtraction

        strand = other.strand_valid if strand_behavior == STRAND_BEHAVIOR_AUTO else strand_behavior != "ignore"

        other_clusters = other.merge_overlaps(use_strand=strand, match_by=match_by)

        _by = group_keys_from_strand_behavior(self, other_clusters, strand_behavior=strand_behavior, by=match_by)

        gr = self.count_overlaps(
            other_clusters,
            strand_behavior=strand_behavior,
            overlap_col=TEMP_NUM_COL,
            match_by=_by,
        )

        result = gr.apply_pair(
            other_clusters,
            strand_behavior=strand_behavior,
            function=_subtraction,
            by=_by,
            skip_if_empty=False,
        )

        return result.drop_and_return(TEMP_NUM_COL, axis=1)

    def summary(
        self,
        *,
        return_df: bool = False,
    ) -> pd.DataFrame | None:
        """Return info.

        Count refers to the number of intervals, the rest to the lengths.

        The column "pyrange" describes the data as is. "coverage_forward" and "coverage_reverse"
        describe the data after strand-specific merging of overlapping intervals.
        "coverage_unstranded" describes the data after merging, without considering the strands.

        The row "count" is the number of intervals and "sum" is their total length. The rest describe the lengths of the
        intervals.

        Parameters
        ----------
        return_df : bool, default False
            Return df with summary.

        Returns
        -------
            None or pd.DataFrame with summary.


        Examples
        --------
        >>> gr = pr.example_data.ensembl_gtf.get_with_loc_columns(["Feature", "gene_id"])
        >>> gr
        index    |    Chromosome    Start    End      Strand      Feature     gene_id
        int64    |    category      int64    int64    category    category    object
        -------  ---  ------------  -------  -------  ----------  ----------  ---------------
        0        |    1             11868    14409    +           gene        ENSG00000223972
        1        |    1             11868    14409    +           transcript  ENSG00000223972
        2        |    1             11868    12227    +           exon        ENSG00000223972
        3        |    1             12612    12721    +           exon        ENSG00000223972
        ...      |    ...           ...      ...      ...         ...         ...
        7        |    1             120724   133723   -           transcript  ENSG00000238009
        8        |    1             133373   133723   -           exon        ENSG00000238009
        9        |    1             129054   129223   -           exon        ENSG00000238009
        10       |    1             120873   120932   -           exon        ENSG00000238009
        PyRanges with 11 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.summary()
                 pyrange    coverage_forward    coverage_reverse    coverage_unstranded
        -----  ---------  ------------------  ------------------  ---------------------
        count      11                      1                3                      4
        mean     1893.27                2541             4503                   4012.5
        std      3799.24                 nan             7359.28                6088.38
        min        59                   2541              105                    105
        25%       139                   2541              255                    330
        50%       359                   2541              405                   1473
        75%      1865                   2541             6702                   5155.5
        max     12999                   2541            12999                  12999
        sum     20826                   2541            13509                  16050

        >>> gr.summary(return_df=True)
                    pyrange  coverage_forward  coverage_reverse  coverage_unstranded
        count     11.000000               1.0          3.000000             4.000000
        mean    1893.272727            2541.0       4503.000000          4012.500000
        std     3799.238610               NaN       7359.280671          6088.379834
        min       59.000000            2541.0        105.000000           105.000000
        25%      139.000000            2541.0        255.000000           330.000000
        50%      359.000000            2541.0        405.000000          1473.000000
        75%     1865.000000            2541.0       6702.000000          5155.500000
        max    12999.000000            2541.0      12999.000000         12999.000000
        sum    20826.000000            2541.0      13509.000000         16050.000000

        """
        from pyranges.methods.summary import _summary

        return _summary(self, return_df=return_df)

    def tile(
        self,
        tile_size: int,
        *,
        overlap_column: str | None = None,
    ) -> "PyRanges":
        """Return overlapping genomic tiles.

        The genome is divided into bookended tiles of length `tile_size`. One tile is returned for each
        interval that overlaps with it, including any metadata from the original intervals.

        Parameters
        ----------
        tile_size : int
            Length of the tiles.

        overlap_column
            Name of column to add with the overlap between each bookended tile.



        Returns
        -------
        PyRanges

            Tiled PyRanges.

        Warning
        -------
        The returned Pyranges may have index duplicates. Call .reset_index(drop=True) to fix it.



        See Also
        --------
        PyRanges.window : divide intervals into windows

        Examples
        --------
        >>> gr = pr.example_data.ensembl_gtf.get_with_loc_columns(["Feature", "gene_name"])
        >>> gr
        index    |    Chromosome    Start    End      Strand      Feature     gene_name
        int64    |    category      int64    int64    category    category    object
        -------  ---  ------------  -------  -------  ----------  ----------  -----------
        0        |    1             11868    14409    +           gene        DDX11L1
        1        |    1             11868    14409    +           transcript  DDX11L1
        2        |    1             11868    12227    +           exon        DDX11L1
        3        |    1             12612    12721    +           exon        DDX11L1
        ...      |    ...           ...      ...      ...         ...         ...
        7        |    1             120724   133723   -           transcript  AL627309.1
        8        |    1             133373   133723   -           exon        AL627309.1
        9        |    1             129054   129223   -           exon        AL627309.1
        10       |    1             120873   120932   -           exon        AL627309.1
        PyRanges with 11 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.tile(200)
        index    |    Chromosome    Start    End      Strand      Feature     gene_name
        int64    |    category      int64    int64    category    category    object
        -------  ---  ------------  -------  -------  ----------  ----------  -----------
        0        |    1             11800    12000    +           gene        DDX11L1
        0        |    1             12000    12200    +           gene        DDX11L1
        0        |    1             12200    12400    +           gene        DDX11L1
        0        |    1             12400    12600    +           gene        DDX11L1
        ...      |    ...           ...      ...      ...         ...         ...
        8        |    1             133600   133800   -           exon        AL627309.1
        9        |    1             129000   129200   -           exon        AL627309.1
        9        |    1             129200   129400   -           exon        AL627309.1
        10       |    1             120800   121000   -           exon        AL627309.1
        PyRanges with 116 rows, 6 columns, and 1 index columns (with 105 index duplicates).
        Contains 1 chromosomes and 2 strands.

        >>> gr.tile(100, overlap_column="TileOverlap")
        index    |    Chromosome    Start    End      Strand      Feature     gene_name    TileOverlap
        int64    |    category      int64    int64    category    category    object       int64
        -------  ---  ------------  -------  -------  ----------  ----------  -----------  -------------
        0        |    1             11800    11900    +           gene        DDX11L1      32
        0        |    1             11900    12000    +           gene        DDX11L1      100
        0        |    1             12000    12100    +           gene        DDX11L1      100
        0        |    1             12100    12200    +           gene        DDX11L1      100
        ...      |    ...           ...      ...      ...         ...         ...          ...
        9        |    1             129100   129200   -           exon        AL627309.1   100
        9        |    1             129200   129300   -           exon        AL627309.1   23
        10       |    1             120800   120900   -           exon        AL627309.1   27
        10       |    1             120900   121000   -           exon        AL627309.1   32
        PyRanges with 223 rows, 7 columns, and 1 index columns (with 212 index duplicates).
        Contains 1 chromosomes and 2 strands.

        >>> gr.tile(100, overlap_column="TileOverlap").reset_index(drop=True)
        index    |    Chromosome    Start    End      Strand      Feature     gene_name    TileOverlap
        int64    |    category      int64    int64    category    category    object       int64
        -------  ---  ------------  -------  -------  ----------  ----------  -----------  -------------
        0        |    1             11800    11900    +           gene        DDX11L1      32
        1        |    1             11900    12000    +           gene        DDX11L1      100
        2        |    1             12000    12100    +           gene        DDX11L1      100
        3        |    1             12100    12200    +           gene        DDX11L1      100
        ...      |    ...           ...      ...      ...         ...         ...          ...
        219      |    1             129100   129200   -           exon        AL627309.1   100
        220      |    1             129200   129300   -           exon        AL627309.1   23
        221      |    1             120800   120900   -           exon        AL627309.1   27
        222      |    1             120900   121000   -           exon        AL627309.1   32
        PyRanges with 223 rows, 7 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.


        """
        from pyranges.methods.windows import _tiles

        kwargs = {
            "overlap_column": overlap_column,
            "tile_size": tile_size,
        }

        # every interval can be processed individually. This may be optimized in the future
        res = self.apply_single(_tiles, by=None, preserve_index=True, **kwargs)
        return mypy_ensure_pyranges(res)

    def three_end(
        self,
        transcript_id: str | list[str] | None = None,
    ) -> "PyRanges":
        """Return the 3'-end.

        The 3'-end is the the end of intervals on the forward strand and the start of intervals on the reverse strand.

        Parameters
        ----------
        transcript_id : str or list of str, default: None
            Optional column name(s). If provided, the three prime end is calculated for each
            group of intervals.


        Returns
        -------
        PyRanges
            PyRanges with the three prime ends

        See Also
        --------
        PyRanges.five_end : return the five prime end
        PyRanges.subsequence : return subintervals specified in relative genome-based coordinates
        PyRanges.spliced_subsequence : return subintervals specified in relative mRNA-based coordinates

        Examples
        --------
        >>> gr = pr.PyRanges({'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 10, 5], 'End': [9, 14, 7],
        ...                    'Strand': ["+", "+", "-"], 'Name': ['a', 'a', 'b']})
        >>> gr
          index  |    Chromosome      Start      End  Strand    Name
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  --------
              0  |    chr1                3        9  +         a
              1  |    chr1               10       14  +         a
              2  |    chr1                5        7  -         b
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.three_end()
          index  |    Chromosome      Start      End  Strand    Name
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  --------
              0  |    chr1                8        9  +         a
              1  |    chr1               13       14  +         a
              2  |    chr1                5        6  -         b
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.three_end(transcript_id='Name')
          index  |    Chromosome      Start      End  Strand    Name
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  --------
              1  |    chr1               13       14  +         a
              2  |    chr1                6        7  -         b
        PyRanges with 2 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        """
        if not (self.has_strand and self.strand_valid):
            msg = f"Need PyRanges with valid strands ({VALID_GENOMIC_STRAND_INFO}) to find 3'."
            raise AssertionError(msg)
        return (
            mypy_ensure_pyranges(
                self.apply_single(function=_tes, by=None, use_strand=True),
            )
            if transcript_id is None
            else self.subsequence(transcript_id=transcript_id, start=-1)
        )

    def to_bed(
        self,
        path: str | None = None,
        compression: PANDAS_COMPRESSION_TYPE = None,
        *,
        keep: bool = True,
    ) -> str | None:
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
          index  |    Chromosome      Start      End  Strand       Gene
          int64  |    object          int64    int64  object      int64
        -------  ---  ------------  -------  -------  --------  -------
              0  |    chr1                1        5  +               1
              1  |    chr1                6        8  -               2
        PyRanges with 2 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

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
        from pyranges.core.out import _to_bed

        return _to_bed(self, path, keep=keep, compression=compression)

    def to_bigwig(
        self: "pr.PyRanges",
        path: None = None,
        chromosome_sizes: pd.DataFrame | dict | None = None,
        value_col: str | None = None,
        *,
        divide: bool = False,
        rpm: bool = True,
        dryrun: bool = False,
        chain: bool = False,
    ) -> "PyRanges | None":
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

        chain: bool, default False
            Return the bigwig data created

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
          index  |    Chromosome      Start      End  Strand      Value
          int64  |    object          int64    int64  object      int64
        -------  ---  ------------  -------  -------  --------  -------
              0  |    chr1                1        7  +              10
              1  |    chr1                4        8  -              20
              2  |    chr1                6       10  -              30
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.to_bigwig(dryrun=True, rpm=False)
          index  |    Chromosome      Start      End      Score
          int64  |    category        int64    int64    float64
        -------  ---  ------------  -------  -------  ---------
              1  |    chr1                1        4          1
              2  |    chr1                4        6          2
              3  |    chr1                6        7          3
              4  |    chr1                7        8          2
              5  |    chr1                8       10          1
        PyRanges with 5 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.to_bigwig(dryrun=True, rpm=False, value_col="Value")
          index  |    Chromosome      Start      End      Score
          int64  |    category        int64    int64    float64
        -------  ---  ------------  -------  -------  ---------
              1  |    chr1                1        4         10
              2  |    chr1                4        6         30
              3  |    chr1                6        7         60
              4  |    chr1                7        8         50
              5  |    chr1                8       10         30
        PyRanges with 5 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.to_bigwig(dryrun=True, rpm=False, value_col="Value", divide=True)
          index  |    Chromosome      Start      End      Score
          int64  |    category        int64    int64    float64
        -------  ---  ------------  -------  -------  ---------
              0  |    chr1                0        1  nan
              1  |    chr1                1        4    3.32193
              2  |    chr1                4        6    3.90689
              3  |    chr1                6        7    4.32193
              4  |    chr1                7        8    4.64386
              5  |    chr1                8       10    4.90689
        PyRanges with 6 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        from pyranges.core.out import _to_bigwig

        _chromosome_sizes = pr.example_data.chromsizes if chromosome_sizes is None else chromosome_sizes

        result = _to_bigwig(
            self=self,
            path=path,
            chromosome_sizes=_chromosome_sizes,
            rpm=rpm,
            divide=divide,
            value_col=value_col,
            dryrun=dryrun,
        )

        if dryrun:
            return result

        if chain:
            return self
        return None

    def to_gff3(
        self,
        path: None = None,
        compression: PANDAS_COMPRESSION_TYPE = None,
    ) -> str | None:
        r"""Write to General Feature Format 3.

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

        compression : {'infer', 'gzip', 'bz2', 'zip', 'xz', None}, default "infer"
            Which compression to use. Uses file extension to infer by default.

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
        '1\t.\tgene\t2\t4\t.\t.\t.\t\n1\t.\texon\t4\t6\t.\t.\t.\t\n1\t.\texon\t6\t9\t.\t.\t.\t\n'

        # How the file would look
        1	.	gene	2	4	.	.	.
        1	.	exon	4	6	.	.	.
        1	.	exon	6	9	.	.	.

        >>> gr["Gene"] = [1, 2, 3]
        >>> gr["function"] = ["a b", "c", "def"]
        >>> gr.to_gff3()
        '1\t.\tgene\t2\t4\t.\t.\t.\tGene=1;function=a b\n1\t.\texon\t4\t6\t.\t.\t.\tGene=2;function=c\n1\t.\texon\t6\t9\t.\t.\t.\tGene=3;function=def\n'

        # How the file would look
        1	.	gene	2	4	.	.	.	Gene=1;function=a b
        1	.	exon	4	6	.	.	.	Gene=2;function=c
        1	.	exon	6	9	.	.	.	Gene=3;function=def

        >>> gr["phase"] = [0, 2, 1]
        >>> gr["Feature"] = ['mRNA', 'CDS', 'CDS']
        >>> gr
          index  |      Chromosome    Start      End  Feature       Gene  function      phase
          int64  |           int64    int64    int64  object       int64  object        int64
        -------  ---  ------------  -------  -------  ---------  -------  ----------  -------
              0  |               1        1        4  mRNA             1  a b               0
              1  |               1        3        6  CDS              2  c                 2
              2  |               1        5        9  CDS              3  def               1
        PyRanges with 3 rows, 7 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.to_gff3()
        '1\t.\tmRNA\t2\t4\t.\t.\t0\tGene=1;function=a b\n1\t.\tCDS\t4\t6\t.\t.\t2\tGene=2;function=c\n1\t.\tCDS\t6\t9\t.\t.\t1\tGene=3;function=def\n'

        # How the file would look
        1	.	mRNA	2	4	.	.	0	Gene=1;function=a b
        1	.	CDS	4	6	.	.	2	Gene=2;function=c
        1	.	CDS	6	9	.	.	1	Gene=3;function=def

        """
        from pyranges.core.out import _to_gff_like

        return _to_gff_like(self, path=path, out_format="gff3", compression=compression)

    def to_gtf(
        self,
        path: None = None,
        compression: PANDAS_COMPRESSION_TYPE = None,
    ) -> str | None:
        r"""Write to Gene Transfer Format.

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

        compression : {'infer', 'gzip', 'bz2', 'zip', 'xz', None}, default "infer"
            Which compression to use. Uses file extension to infer by default.

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
        '1\t.\tgene\t2\t4\t.\t.\t.\t\n1\t.\texon\t4\t6\t.\t.\t.\t\n1\t.\texon\t6\t9\t.\t.\t.\t\n'

        # What the file contents look like:
        1	.	gene	2	4	.	.	.
        1	.	exon	4	6	.	.	.
        1	.	exon	6	9	.	.	.

        >>> gr.Feature = ["GENE", "EXON", "EXON"]
        >>> gr.to_gtf()  # the raw string output
        '1\t.\tGENE\t2\t4\t.\t.\t.\t\n1\t.\tEXON\t4\t6\t.\t.\t.\t\n1\t.\tEXON\t6\t9\t.\t.\t.\t\n'

        The file would look like:

            1	.	GENE	2	4	.	.	.
            1	.	EXON	4	6	.	.	.
            1	.	EXON	6	9	.	.	.

        """
        from pyranges.core.out import _to_gff_like

        return _to_gff_like(self, path=path, out_format="gtf", compression=compression)

    def to_rle(
        self,
        value_col: str | None = None,
        strand: VALID_USE_STRAND_TYPE = "auto",
        *,
        rpm: bool = False,
    ) -> "Rledict":
        """Return as Rledict.

        Create collection of Rles representing the coverage or other numerical value.

        Parameters
        ----------
        value_col : str, default None
            Numerical column to create Rledict from.

        strand : bool, default None, i.e. auto
            Whether to treat strands serparately.

        rpm : bool, default False
            Normalize by multiplying with `1e6/(number_intervals)`.

        Returns
        -------
        pyrle.Rledict

            Rle with coverage or other info from the PyRanges.

        Examples
        --------
        # >>> d = {'Chromosome': ['chr1', 'chr1', 'chr1'], 'Start': [3, 8, 5],
        # ...      'End': [6, 9, 7], 'Score': [0.1, 5, 3.14], 'Strand': ['+', '+', '-']}
        # >>> gr = pr.PyRanges(d)
        # >>> gr.to_rle()
        # chr1 +
        # --
        # +--------+-----+-----+-----+-----+
        # | Runs   | 3   | 3   | 2   | 1   |
        # |--------+-----+-----+-----+-----|
        # | Values | 0.0 | 1.0 | 0.0 | 1.0 |
        # +--------+-----+-----+-----+-----+
        # Rle of length 9 containing 4 elements (avg. length 2.25)
        # <BLANKLINE>
        # chr1 -
        # --
        # +--------+-----+-----+
        # | Runs   | 5   | 2   |
        # |--------+-----+-----|
        # | Values | 0.0 | 1.0 |
        # +--------+-----+-----+
        # Rle of length 7 containing 2 elements (avg. length 3.5)
        # RleDict object with 2 chromosomes/strand pairs.

        # >>> gr.to_rle(value_col="Score")
        # chr1 +
        # --
        # +--------+-----+-----+-----+-----+
        # | Runs   | 3   | 3   | 2   | 1   |
        # |--------+-----+-----+-----+-----|
        # | Values | 0.0 | 0.1 | 0.0 | 5.0 |
        # +--------+-----+-----+-----+-----+
        # Rle of length 9 containing 4 elements (avg. length 2.25)
        # <BLANKLINE>
        # chr1 -
        # --
        # +--------+-----+------+
        # | Runs   | 5   | 2    |
        # |--------+-----+------|
        # | Values | 0.0 | 3.14 |
        # +--------+-----+------+
        # Rle of length 7 containing 2 elements (avg. length 3.5)
        # RleDict object with 2 chromosomes/strand pairs.

        # >>> gr.to_rle(value_col="Score", strand=False)
        # chr1
        # +--------+-----+-----+------+------+-----+-----+
        # | Runs   | 3   | 2   | 1    | 1    | 1   | 1   |
        # |--------+-----+-----+------+------+-----+-----|
        # | Values | 0.0 | 0.1 | 3.24 | 3.14 | 0.0 | 5.0 |
        # +--------+-----+-----+------+------+-----+-----+
        # Rle of length 9 containing 6 elements (avg. length 1.5)
        # Unstranded Rledict object with 1 chromosome.

        # >>> gr.to_rle(rpm=True)
        # chr1 +
        # --
        # +--------+-----+-------------------+-----+-------------------+
        # | Runs   | 3   | 3                 | 2   | 1                 |
        # |--------+-----+-------------------+-----+-------------------|
        # | Values | 0.0 | 333333.3333333333 | 0.0 | 333333.3333333333 |
        # +--------+-----+-------------------+-----+-------------------+
        # Rle of length 9 containing 4 elements (avg. length 2.25)
        # <BLANKLINE>
        # chr1 -
        # --
        # +--------+-----+-------------------+
        # | Runs   | 5   | 2                 |
        # |--------+-----+-------------------|
        # | Values | 0.0 | 333333.3333333333 |
        # +--------+-----+-------------------+
        # Rle of length 7 containing 2 elements (avg. length 3.5)
        # Rledict object with 2 chromosomes/strand pairs.

        """
        if strand is None:
            strand = self.strand_valid

        from pyranges.methods.to_rle import _to_rle

        strand = validate_and_convert_strand(self, strand)
        df = self.remove_strand() if not strand else self

        return _to_rle(
            df,
            value_col,
            strand=strand,
            rpm=rpm,
        )

    def remove_strand(self) -> "PyRanges":
        """Return a copy with the Strand column removed.

        Strand is removed regardless of whether it contains valid strand info.

        See Also
        --------
        PyRanges.strand_valid : whether PyRanges has valid strand info
        PyRanges.remove_strand : remove the Strand column from PyRanges

        Examples
        --------
        >>> gr = pr.PyRanges({'Chromosome': ['chr1', 'chr1'], 'Start': [1, 6],
        ...                   'End': [5, 8], 'Strand': ['+', '-']})
        >>> gr
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        5  +
              1  |    chr1                6        8  -
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.remove_strand()
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                1        5
              1  |    chr1                6        8
        PyRanges with 2 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        if not self.has_strand:
            return self
        return self.drop_and_return(STRAND_COL, axis=1)

    def window(
        self,
        window_size: int,
        use_strand: VALID_USE_STRAND_TYPE = USE_STRAND_DEFAULT,
    ) -> "PyRanges":
        """Return non-overlapping genomic windows.

        Every interval is split into windows of length `window_size` starting from its 5' end.

        Parameters
        ----------
        window_size : int
            Length of the windows.

        use_strand: {"auto", True, False}, default: "auto"
            Whether negative strand intervals should be sliced in descending order, meaning 5' to 3'.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        Returns
        -------
        PyRanges

            Sliding window PyRanges.

        Warning
        -------
        The returned Pyranges may have index duplicates. Call .reset_index(drop=True) to fix it.

        See Also
        --------
        PyRanges.tile : divide intervals into adjacent tiles.

        Examples
        --------
        >>> import pyranges as pr
        >>> gr = pr.PyRanges({"Chromosome": [1], "Start": [800], "End": [1012]})
        >>> gr
          index  |      Chromosome    Start      End
          int64  |           int64    int64    int64
        -------  ---  ------------  -------  -------
              0  |               1      800     1012
        PyRanges with 1 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.window(100)
          index  |      Chromosome    Start      End
          int64  |           int64    int64    int64
        -------  ---  ------------  -------  -------
              0  |               1      800      900
              0  |               1      900     1000
              0  |               1     1000     1012
        PyRanges with 3 rows, 3 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes.

        >>> gr.window(100).reset_index(drop=True)
          index  |      Chromosome    Start      End
          int64  |           int64    int64    int64
        -------  ---  ------------  -------  -------
              0  |               1      800      900
              1  |               1      900     1000
              2  |               1     1000     1012
        PyRanges with 3 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        # Negative strand intervals are sliced in descending order by default:
        >>> gs = pr.PyRanges({"Chromosome": [1, 1], "Start": [200, 600], "End": [332, 787], "Strand":['+', '-']})
        >>> gs
          index  |      Chromosome    Start      End  Strand
          int64  |           int64    int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |               1      200      332  +
              1  |               1      600      787  -
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> w=gs.window(100)
        >>> w['lengths']=w.lengths() # add lengths column to see the length of the windows
        >>> w
          index  |      Chromosome    Start      End  Strand      lengths
          int64  |           int64    int64    int64  object        int64
        -------  ---  ------------  -------  -------  --------  ---------
              0  |               1      200      300  +               100
              0  |               1      300      332  +                32
              1  |               1      687      787  -               100
              1  |               1      600      687  -                87
        PyRanges with 4 rows, 5 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes and 2 strands.

        >>> gs.window(100, use_strand=False)
          index  |      Chromosome    Start      End  Strand
          int64  |           int64    int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |               1      200      300  +
              0  |               1      300      332  +
              1  |               1      600      700  -
              1  |               1      700      787  -
        PyRanges with 4 rows, 4 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes and 2 strands.

        >>> gr2 = pr.example_data.ensembl_gtf.get_with_loc_columns(["Feature", "gene_name"])
        >>> gr2
        index    |    Chromosome    Start    End      Strand      Feature     gene_name
        int64    |    category      int64    int64    category    category    object
        -------  ---  ------------  -------  -------  ----------  ----------  -----------
        0        |    1             11868    14409    +           gene        DDX11L1
        1        |    1             11868    14409    +           transcript  DDX11L1
        2        |    1             11868    12227    +           exon        DDX11L1
        3        |    1             12612    12721    +           exon        DDX11L1
        ...      |    ...           ...      ...      ...         ...         ...
        7        |    1             120724   133723   -           transcript  AL627309.1
        8        |    1             133373   133723   -           exon        AL627309.1
        9        |    1             129054   129223   -           exon        AL627309.1
        10       |    1             120873   120932   -           exon        AL627309.1
        PyRanges with 11 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr2 = pr.example_data.ensembl_gtf.get_with_loc_columns(["Feature", "gene_name"])
        >>> gr2.window(1000)
        index    |    Chromosome    Start    End      Strand      Feature     gene_name
        int64    |    category      int64    int64    category    category    object
        -------  ---  ------------  -------  -------  ----------  ----------  -----------
        0        |    1             11868    12868    +           gene        DDX11L1
        0        |    1             12868    13868    +           gene        DDX11L1
        0        |    1             13868    14409    +           gene        DDX11L1
        1        |    1             11868    12868    +           transcript  DDX11L1
        ...      |    ...           ...      ...      ...         ...         ...
        7        |    1             120724   121723   -           transcript  AL627309.1
        8        |    1             133373   133723   -           exon        AL627309.1
        9        |    1             129054   129223   -           exon        AL627309.1
        10       |    1             120873   120932   -           exon        AL627309.1
        PyRanges with 28 rows, 6 columns, and 1 index columns (with 17 index duplicates).
        Contains 1 chromosomes and 2 strands.

        """
        from pyranges.methods.windows import _windows

        use_strand = validate_and_convert_strand(self, use_strand)

        kwargs = {
            "window_size": window_size,
        }

        # every interval can be processed individually. This may be optimized in the future.
        df = self.apply_single(_windows, by=None, use_strand=use_strand, preserve_index=True, **kwargs)
        return mypy_ensure_pyranges(df)

    def remove_nonloc_columns(self) -> "PyRanges":
        """Remove all columns that are not genome location columns (Chromosome, Start, End, Strand).

        Examples
        --------
        >>> gr = pr.PyRanges({"Chromosome": [1], "Start": [895], "Strand": ["+"], "Score": [1], "Score2": [2], "End": [1259]})
        >>> gr
          index  |      Chromosome    Start  Strand      Score    Score2      End
          int64  |           int64    int64  object      int64     int64    int64
        -------  ---  ------------  -------  --------  -------  --------  -------
              0  |               1      895  +               1         2     1259
        PyRanges with 1 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.
        >>> gr.remove_nonloc_columns()
          index  |      Chromosome    Start      End  Strand
          int64  |           int64    int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |               1      895     1259  +
        PyRanges with 1 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        """
        cols = GENOME_LOC_COLS_WITH_STRAND if self.has_strand else GENOME_LOC_COLS
        return mypy_ensure_pyranges(super().__getitem__(cols))

    def get_with_loc_columns(
        self,
        key: str | Iterable[str],
        *,
        preserve_loc_order: bool = False,
    ) -> "pr.PyRanges":
        """Return a PyRanges with the requested columns, as well as the genome location columns.

        Parameters
        ----------
        key : str or iterable of str
            Column(s) to return

        preserve_loc_order : bool, default False
            Whether to preserve the order of the genome location columns.
            If False, the genome location columns will be moved to the left.

        Returns
        -------
        PyRanges

            PyRanges with the requested columns.

        >>> gr = pr.PyRanges({"Chromosome": [1], "Start": [895], "Strand": ["+"], "Score": [1], "Score2": [2], "End": [1259]})
        >>> gr
          index  |      Chromosome    Start  Strand      Score    Score2      End
          int64  |           int64    int64  object      int64     int64    int64
        -------  ---  ------------  -------  --------  -------  --------  -------
              0  |               1      895  +               1         2     1259
        PyRanges with 1 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        >>> gr.get_with_loc_columns(["Score2", "Score", "Score2"]) # moves loc columns to the left by default
          index  |      Chromosome    Start      End  Strand      Score2    Score    Score2
          int64  |           int64    int64    int64  object       int64    int64     int64
        -------  ---  ------------  -------  -------  --------  --------  -------  --------
              0  |               1      895     1259  +                2        1         2
        PyRanges with 1 rows, 7 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        >>> gr.get_with_loc_columns(["Score2", "Score"], preserve_loc_order=True)
          index  |      Chromosome    Start  Strand      Score2    Score      End
          int64  |           int64    int64  object       int64    int64    int64
        -------  ---  ------------  -------  --------  --------  -------  -------
              0  |               1      895  +                2        1     1259
        PyRanges with 1 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        >>> gr.get_with_loc_columns(["Score2", "Score", "Score2"], preserve_loc_order=True)
        Traceback (most recent call last):
        ...
        ValueError: Duplicate keys not allowed when preserve_loc_order is True.

        """
        keys = [key] if isinstance(key, str) else ([*key])

        def _reorder_according_to_b(a: list[str], b: list[str]) -> list[str]:
            for pos, val in zip(sorted([a.index(x) for x in b]), b, strict=True):
                a[pos] = val
            return a

        if preserve_loc_order:
            if len(set(keys)) != len(keys):
                msg = "Duplicate keys not allowed when preserve_loc_order is True."
                raise ValueError(msg)
            cols_to_include = {*keys, *GENOME_LOC_COLS_WITH_STRAND}
            cols_to_include_genome_loc_correct_order = [col for col in self.columns if col in cols_to_include]
            cols_to_include_genome_loc_correct_order = _reorder_according_to_b(
                cols_to_include_genome_loc_correct_order,
                keys,
            )
        else:
            loc_columns = GENOME_LOC_COLS_WITH_STRAND if self.has_strand else GENOME_LOC_COLS
            cols_to_include_genome_loc_correct_order = loc_columns + keys

        return mypy_ensure_pyranges(super().__getitem__(cols_to_include_genome_loc_correct_order))

    # TO DO: unclear how argument
    def intersect(  # type: ignore[override]
        self,
        other: "PyRanges",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        multiple: VALID_OVERLAP_TYPE = "all",
        *,
        match_by: VALID_BY_TYPES = None,
        **_,
    ) -> "PyRanges":
        """Return overlapping subintervals.

        Returns the segments of the intervals in self which overlap with those in other.
        When multiple intervals in 'other' overlap with the same interval in self, the result
        may be complex -- read the argument 'multiple' for details.

        Parameters
        ----------
        other : PyRanges
            PyRanges to find overlaps with.

        multiple : {"all", "first", "last"}, default "all"
            What intervals to report when multiple intervals in 'other' overlap with the same interval in self.
            The default "all" reports all overlapping subintervals, which will have duplicate indices.
            "first" reports only, for each interval in self, the overlapping subinterval with smallest Start in 'other'
            "last" reports only the overlapping subinterval with the biggest End in 'other'

        strand_behavior : {"auto", "same", "opposite", "ignore"}, default "auto"
            Whether to consider overlaps of intervals on the same strand, the opposite or ignore strand
            information. The default, "auto", means use "same" if both PyRanges are stranded (see .strand_valid)
            otherwise ignore the strand information.

        match_by : str or list, default None
            If provided, only intervals with an equal value in column(s) `match_by` may be considered as overlapping.

        Returns
        -------
        PyRanges

            A PyRanges with overlapping intervals.

        See Also
        --------
        PyRanges.overlap : report overlapping (unmodified) intervals
        PyRanges.subtract_ranges : report non-overlapping subintervals
        PyRanges.set_intersect : set-intersect PyRanges

        Examples
        --------
        >>> r1 = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [5, 20, 40],"End": [10, 30, 50], "ID": ["a", "b", "c"]})
        >>> r1
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                5       10  a
              1  |    chr1               20       30  b
              2  |    chr1               40       50  c
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.


        >>> r2 = pr.PyRanges({"Chromosome": ["chr1"] * 4, "Start": [7, 18, 25, 28], "End": [9, 22, 33, 32]})
        >>> r2
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                7        9
              1  |    chr1               18       22
              2  |    chr1               25       33
              3  |    chr1               28       32
        PyRanges with 4 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> r1.intersect(r2)
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                7        9  a
              1  |    chr1               20       22  b
              1  |    chr1               25       30  b
              1  |    chr1               28       30  b
        PyRanges with 4 rows, 4 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes.

        >>> r1.intersect(r2, multiple="first")
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                7        9  a
              1  |    chr1               20       22  b
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> r1.intersect(r2, multiple="last")
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                7        9  a
              1  |    chr1               25       30  b
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        from pyranges.methods.overlap import _intersect

        # note: argument multiple = 'containment' is formally accepted but omitted in docstring since the result
        # will be always the same as self.overlap(other, contained=True), no intersect is done in that case

        return self.apply_pair(
            other,
            _intersect,
            strand_behavior=strand_behavior,
            by=match_by,
            how=multiple,
        )

    def combine_interval_columns(
        self,
        function: VALID_COMBINE_OPTIONS | CombineIntervalColumnsOperation = "intersect",
        *,
        start: str = START_COL,
        end: str = END_COL,
        start2: str = START_COL + JOIN_SUFFIX,
        end2: str = END_COL + JOIN_SUFFIX,
        drop_old_columns: bool = True,
    ) -> "pr.PyRanges":
        """Use two pairs of columns representing intervals to create a new start and end column.

        By default, the new start and end columns will be the intersection of the intervals.

        Parameters
        ----------
        function : {"intersect", "union"} or Callable, default "intersect"
            How to combine the intervals: "intersect" or "union".
            If a callable is passed, it should take four Series arguments: start1, end1, start2, end2;
            and return a tuple of two integers: (new_starts, new_ends).

        start : str, default "Start"
            Column name for Start of first interval
        end : str, default "End"
            Column name for End of first interval
        start2 : str, default "Start_b"
            Column name for Start of second interval
        end2 : str, default "End_b"
            Column name for End of second interval
        drop_old_columns : bool, default True
            Whether to drop the above mentioned columns.

        Examples
        --------
        >>> gr1, gr2 = pr.example_data.aorta.head(3).remove_nonloc_columns(), pr.example_data.aorta2.head(3).remove_nonloc_columns()
        >>> j = gr1.join_ranges(gr2)
        >>> j
          index  |    Chromosome      Start      End  Strand      Chromosome_b      Start_b    End_b  Strand_b
          int64  |    category        int64    int64  category    category            int64    int64  category
        -------  ---  ------------  -------  -------  ----------  --------------  ---------  -------  ----------
              1  |    chr1             9939    10138  +           chr1                10073    10272  +
              0  |    chr1             9916    10115  -           chr1                 9988    10187  -
              0  |    chr1             9916    10115  -           chr1                10079    10278  -
              2  |    chr1             9951    10150  -           chr1                 9988    10187  -
              2  |    chr1             9951    10150  -           chr1                10079    10278  -
        PyRanges with 5 rows, 8 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes and 2 strands.

        # intersect the intervals by default
        >>> j.combine_interval_columns()
          index  |    Chromosome      Start      End  Strand      Chromosome_b    Strand_b
          int64  |    category        int64    int64  category    category        category
        -------  ---  ------------  -------  -------  ----------  --------------  ----------
              1  |    chr1            10073    10138  +           chr1            +
              0  |    chr1             9988    10115  -           chr1            -
              0  |    chr1            10079    10115  -           chr1            -
              2  |    chr1             9988    10150  -           chr1            -
              2  |    chr1            10079    10150  -           chr1            -
        PyRanges with 5 rows, 6 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes and 2 strands.

        # take the union instead
        >>> j.combine_interval_columns('union')
          index  |    Chromosome      Start      End  Strand      Chromosome_b    Strand_b
          int64  |    category        int64    int64  category    category        category
        -------  ---  ------------  -------  -------  ----------  --------------  ----------
              1  |    chr1             9939    10272  +           chr1            +
              0  |    chr1             9916    10187  -           chr1            -
              0  |    chr1             9916    10278  -           chr1            -
              2  |    chr1             9951    10187  -           chr1            -
              2  |    chr1             9951    10278  -           chr1            -
        PyRanges with 5 rows, 6 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes and 2 strands.

        # use a custom function that keeps the start of the first interval and the end of the second
        >>> def custom_combine(s1, e1, s2, e2): return (s1, e2)
        >>> j.combine_interval_columns(custom_combine)
          index  |    Chromosome      Start      End  Strand      Chromosome_b    Strand_b
          int64  |    category        int64    int64  category    category        category
        -------  ---  ------------  -------  -------  ----------  --------------  ----------
              1  |    chr1             9939    10272  +           chr1            +
              0  |    chr1             9916    10187  -           chr1            -
              0  |    chr1             9916    10278  -           chr1            -
              2  |    chr1             9951    10187  -           chr1            -
              2  |    chr1             9951    10278  -           chr1            -
        PyRanges with 5 rows, 6 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes and 2 strands.


        """
        from pyranges.methods.combine_positions import _intersect_interval_columns, _union_interval_columns

        if function == "intersect":
            function = _intersect_interval_columns
        elif function == "union":
            function = _union_interval_columns

        new_starts, new_ends = function(self[start], self[end], self[start2], self[end2])

        z = self.copy()
        z[START_COL] = new_starts
        z[END_COL] = new_ends

        cols_to_drop = list({start, end, start2, end2}.difference(RANGE_COLS) if drop_old_columns else {})

        return z.drop_and_return(labels=cols_to_drop, axis="columns")

    def get_sequence(
        self: "PyRanges",
        path: Path | None = None,
        pyfaidx_fasta: Optional["pyfaidx.Fasta"] = None,
    ) -> pd.Series:
        r"""Get the sequence of the intervals from a fasta file.

        Parameters
        ----------
        path : Path
            Path to fasta file. It will be indexed using pyfaidx if an index is not found

        pyfaidx_fasta : pyfaidx.Fasta
            Alternative method to provide fasta target, as a pyfaidx.Fasta object


        Returns
        -------
        Series

            Sequences, one per interval. The series is named 'Sequence'

        Note
        ----

        This function requires the library pyfaidx, it can be installed with
        ``conda install -c bioconda pyfaidx`` or ``pip install pyfaidx``.

        Sorting the PyRanges is likely to improve the speed.
        Intervals on the negative strand will be reverse complemented.

        Warning
        -------

        Note that the names in the fasta header and self.Chromosome must be the same.

        See Also
        --------
        PyRanges.get_transcript_sequence : obtain mRNA sequences, by joining exons belonging to the same transcript


        Examples
        --------
        >>> import pyranges as pr
        >>> gr = pr.PyRanges({"Chromosome": ["chr1", "chr1"],
        ...                   "Start": [5, 0], "End": [8, 5],
        ...                   "Strand": ["+", "-"]})

        >>> gr
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                5        8  +
              1  |    chr1                0        5  -
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> tmp_handle = open("temp.fasta", "w+")
        >>> _ = tmp_handle.write(">chr1\n")
        >>> _ = tmp_handle.write("GTAATCAT\n")
        >>> tmp_handle.close()

        >>> seq = gr.get_sequence("temp.fasta")

        >>> seq
        0      CAT
        1    ATTAC
        Name: Sequence, dtype: object

        >>> gr["seq"] = seq
        >>> gr
          index  |    Chromosome      Start      End  Strand    seq
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  --------
              0  |    chr1                5        8  +         CAT
              1  |    chr1                0        5  -         ATTAC
        PyRanges with 2 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        """
        try:
            import pyfaidx  # type: ignore[import]
        except ImportError:
            LOGGER.exception(
                "pyfaidx must be installed to get fasta sequences. Use `conda install -c bioconda pyfaidx` or `pip install pyfaidx` to install it.",
            )
            sys.exit(1)

        if pyfaidx_fasta is None:
            if path is None:
                msg = "ERROR get_sequence : you must provide a fasta path or pyfaidx_fasta object"
                raise ValueError(msg)
            pyfaidx_fasta = pyfaidx.Fasta(path, read_ahead=int(1e5))

        use_strand = self.strand_valid
        iterables = (
            zip(self[CHROM_COL], self[START_COL], self[END_COL], [FORWARD_STRAND], strict=False)
            if not use_strand
            else zip(self[CHROM_COL], self[START_COL], self[END_COL], self[STRAND_COL], strict=True)
        )
        seqs = []
        for chromosome, start, end, strand in iterables:
            _fasta = pyfaidx_fasta[chromosome]
            forward_strand = strand == FORWARD_STRAND
            if (seq := _fasta[start:end]) is not None:
                seqs.append(seq.seq if forward_strand else (-seq).seq)
        return pd.Series(data=seqs, index=self.index, name="Sequence")

    def get_transcript_sequence(
        self: "PyRanges",
        transcript_id: str,
        path: Path | None = None,
        pyfaidx_fasta: Optional["pyfaidx.Fasta"] = None,
    ) -> pd.DataFrame:
        r"""Get the sequence of mRNAs, e.g. joining intervals corresponding to exons of the same transcript.

        Parameters
        ----------
        transcript_id : str or list of str
            intervals are grouped by this/these ID column(s): these are exons belonging to same transcript

        path : Optional Path
            Path to fasta file. It will be indexed using pyfaidx if an index is not found

        pyfaidx_fasta : pyfaidx.Fasta
            Alternative method to provide fasta target, as a pyfaidx.Fasta object


        Returns
        -------
        DataFrame
            Pandas DataFrame with a column for Sequence, plus ID column(s) provided with "transcript_id"

        Note
        ----

        This function requires the library pyfaidx, it can be installed with
        ``conda install -c bioconda pyfaidx`` or ``pip install pyfaidx``.

        Sorting the PyRanges is likely to improve the speed.
        Intervals on the negative strand will be reverse complemented.

        Warning
        -------

        Note that the names in the fasta header and self.Chromosome must be the same.

        See Also
        --------
        PyRanges.get_sequence : obtain sequence of single intervals


        Examples
        --------
        >>> import pyranges as pr
        >>> gr = pr.PyRanges({"Chromosome": ['chr1'] * 5,
        ...                   "Start": [0, 9, 18, 9, 18], "End": [4, 13, 21, 13, 21],
        ...                   "Strand":['+', '-', '-', '-', '-'],
        ...                   "transcript": ['t1', 't2', 't2', 't4', 't5']})

        >>> tmp_handle = open("temp.fasta", "w+")
        >>> _ = tmp_handle.write(">chr1\n")
        >>> _ = tmp_handle.write("AAACCCTTTGGGAAACCCTTTGGG\n")
        >>> tmp_handle.close()

        >>> seq = gr.get_transcript_sequence(path="temp.fasta", transcript_id='transcript')
        >>> seq
          transcript Sequence
        0         t1     AAAC
        1         t2  AAATCCC
        2         t4     TCCC
        3         t5      AAA

        To write to a file in fasta format:
        >>> with open('outfile.fasta', 'w') as fw:
        ...     nchars=60
        ...     for row in seq.itertuples():
        ...         s = '\\n'.join([ row.Sequence[i:i+nchars] for i in range(0, len(row.Sequence), nchars)])
        ...         _bytes_written = fw.write(f'>{row.transcript}\\n{s}\\n')

        """
        gr = self.sort_by_5_prime_ascending_and_3_prime_descending() if self.strand_valid else self.sort_by_position()

        seq = gr.get_sequence(path=path, pyfaidx_fasta=pyfaidx_fasta)
        gr["Sequence"] = seq.to_numpy()

        return gr.groupby(transcript_id, as_index=False).agg({"Sequence": "".join})

    def genome_bounds(
        self: "PyRanges",
        chromsizes: "dict[str | int, int] | PyRanges",
        *,
        clip: bool = False,
        only_right: bool = False,
    ) -> "PyRanges":
        """Remove or clip intervals outside of genome bounds.

        Parameters
        ----------
        chromsizes : dict or PyRanges or pyfaidx.Fasta
            Dict or PyRanges describing the lengths of the chromosomes.
            pyfaidx.Fasta object is also accepted since it conveniently loads chromosome length

        clip : bool, default False
            Returns the portions of intervals within bounds,
            instead of dropping intervals entirely if they are even partially
            out of bounds

        only_right : bool, default False
            If True, remove or clip only intervals that are out-of-bounds on the right,
            and do not alter those out-of-bounds on the left (whose Start is < 0)


        Examples
        --------
        >>> import pyranges as pr
        >>> d = {"Chromosome": [1, 1, 3], "Start": [1, 249250600, 5], "End": [2, 249250640, 7]}
        >>> gr = pr.PyRanges(d)
        >>> gr
          index  |      Chromosome      Start        End
          int64  |           int64      int64      int64
        -------  ---  ------------  ---------  ---------
              0  |               1          1          2
              1  |               1  249250600  249250640
              2  |               3          5          7
        PyRanges with 3 rows, 3 columns, and 1 index columns.
        Contains 2 chromosomes.

        >>> chromsizes = {1: 249250621, 3: 500}
        >>> chromsizes
        {1: 249250621, 3: 500}

        >>> gr.genome_bounds(chromsizes)
          index  |      Chromosome    Start      End
          int64  |           int64    int64    int64
        -------  ---  ------------  -------  -------
              0  |               1        1        2
              2  |               3        5        7
        PyRanges with 2 rows, 3 columns, and 1 index columns.
        Contains 2 chromosomes.

        >>> gr.genome_bounds(chromsizes, clip=True)
          index  |      Chromosome      Start        End
          int64  |           int64      int64      int64
        -------  ---  ------------  ---------  ---------
              0  |               1          1          2
              1  |               1  249250600  249250621
              2  |               3          5          7
        PyRanges with 3 rows, 3 columns, and 1 index columns.
        Contains 2 chromosomes.

        >>> del chromsizes[3]
        >>> chromsizes
        {1: 249250621}

        >>> gr.genome_bounds(chromsizes)
        Traceback (most recent call last):
        ...
        ValueError: Not all chromosomes were in the chromsize dict. This might mean that their types differed.
        Missing keys: {3}.
        Chromosome col had type: int64 while keys were of type: int

        """
        from pyranges.methods.boundaries import _outside_bounds

        if isinstance(chromsizes, pd.DataFrame):
            chromsizes = dict(*zip(chromsizes[CHROM_COL], chromsizes[END_COL], strict=True))
        elif isinstance(chromsizes, dict):
            pass
        else:  # A hack because pyfaidx might not be installed, but we want type checking anyway
            pyfaidx_chromsizes = cast(dict[str | int, list], chromsizes)
            chromsizes = {k: len(pyfaidx_chromsizes[k]) for k in pyfaidx_chromsizes.keys()}  # noqa: SIM118

        if missing_keys := set(self[CHROM_COL]).difference(set(chromsizes.keys())):
            msg = f"""Not all chromosomes were in the chromsize dict. This might mean that their types differed.
Missing keys: {missing_keys}.
Chromosome col had type: {self[CHROM_COL].dtype} while keys were of type: {', '.join({type(k).__name__ for k in chromsizes})}"""
            raise ValueError(msg)

        if not isinstance(chromsizes, dict):
            msg = "ERROR chromsizes must be a dictionary, or a PyRanges, or a pyfaidx.Fasta object"
            raise TypeError(msg)

        return self.apply_single(
            _outside_bounds,
            by=None,
            use_strand=False,
            preserve_index=True,
            chromsizes=chromsizes,
            clip=clip,
            only_right=only_right,
        )
