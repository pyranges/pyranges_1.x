"""Data structure for genomic intervals and their annotation."""

import logging
import sys
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import TYPE_CHECKING, Optional, cast

import numpy as np
import pandas as pd
from natsort import natsorted  # type: ignore[import]

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
    NEAREST_ANY_DIRECTION,
    NEAREST_DOWNSTREAM,
    NEAREST_UPSTREAM,
    PANDAS_COMPRESSION_TYPE,
    REVERSE_STRAND,
    START_COL,
    STRAND_BEHAVIOR_OPPOSITE,
    STRAND_COL,
    USE_STRAND_DEFAULT,
    VALID_BY_TYPES,
    VALID_COMBINE_OPTIONS,
    VALID_GENOMIC_STRAND_INFO,
    VALID_JOIN_TYPE,
    VALID_NEAREST_TYPE,
    VALID_OVERLAP_TYPE,
    VALID_STRAND_BEHAVIOR_TYPE,
    VALID_USE_STRAND_TYPE,
    CombineIntervalColumnsOperation,
)
from pyranges.core.pyranges_groupby import PyRangesDataFrameGroupBy
from pyranges.core.pyranges_helpers import (
    arg_to_list,
    factorize,
    mypy_ensure_pyranges,
    prepare_by_binary,
    prepare_by_single,
    split_on_strand,
    use_strand_from_validated_strand_behavior,
    validate_and_convert_strand_behavior,
    validate_and_convert_use_strand,
)
from pyranges.core.tostring import tostring
from pyranges.methods.complement import _complement
from pyranges.methods.interval_metrics import compute_interval_metrics
from pyranges.methods.map_to_global import map_to_global
from pyranges.methods.sort import sort_factorize_dict
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

    You can **initialize a PyRanges object like you would a pandas DataFrame**, as long as the resulting DataFrame
    has the necessary columns (Chromosome, Start, End; Strand is optional).
    See examples below, and https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html for more information.

    Parameters
    ----------
    data : dict, pd.DataFrame, or None, default None

    index : Index or array-like

    columns : Index or array-like

    dtyped : type, default None

    copy : bool or None, default None


    See Also
    --------
    pyranges.read_bed: read bed-file into PyRanges
    pyranges.read_bam: read bam-file into PyRanges
    pyranges.read_gff: read gff-file into PyRanges
    pyranges.read_gtf: read gtf-file into PyRanges
    pyranges.from_string: create PyRanges from multiline string

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

    Operations that remove a column required for a PyRanges return a DataFrame instead:

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
        # returning a new instance of a class. It is called before
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

    @property
    def _constructor(self) -> type:
        return pr.PyRanges

    def groupby(self, *args, **kwargs) -> "PyRangesDataFrameGroupBy":
        """Groupby PyRanges."""
        index_of_observed_in_args_list = 7
        if "observed" not in kwargs or len(args) < index_of_observed_in_args_list:
            kwargs["observed"] = True
        grouped = super().groupby(*args, **kwargs)
        return PyRangesDataFrameGroupBy(grouped)

    @property
    def loci(self) -> LociGetter:
        """Get or set rows based on genomic location.

        Parameters
        ----------
        key
            Genomic location: one or more of Chromosome, Strand, and Range (i.e. Start:End).
            When a Range is specified, only rows that overlap with it are returned.

        Returns
        -------
        PyRanges
            PyRanges view with rows matching the location.

        Warning
        ----
        When strand is provided but chromosome is not, only valid strand values ('+', '-') are searched for.
        Use the complete .loci[Chromosome, Strand, Range] syntax to search for non-genomic strands.
        Each item can be None to match all values.

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

        >>> gr.loci["+"]
          index  |      Chromosome    Start      End  Strand      gene_id          gene_name
          int64  |        category    int64    int64  category    object           object
        -------  ---  ------------  -------  -------  ----------  ---------------  -----------
              0  |               1    11868    14409  +           ENSG00000223972  DDX11L1
              1  |               1    11868    14409  +           ENSG00000223972  DDX11L1
              2  |               1    11868    12227  +           ENSG00000223972  DDX11L1
              3  |               1    12612    12721  +           ENSG00000223972  DDX11L1
              4  |               1    13220    14409  +           ENSG00000223972  DDX11L1
        PyRanges with 5 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        >>> gr.loci[11000:12000]
          index  |      Chromosome    Start      End  Strand      gene_id          gene_name
          int64  |        category    int64    int64  category    object           object
        -------  ---  ------------  -------  -------  ----------  ---------------  -----------
              0  |               1    11868    14409  +           ENSG00000223972  DDX11L1
              1  |               1    11868    14409  +           ENSG00000223972  DDX11L1
              2  |               1    11868    12227  +           ENSG00000223972  DDX11L1
        PyRanges with 3 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        The Chromosome column is attempted to be converted to the type of the provided key before matching:

        >>> gr.loci["1", 11000:12000]
          index  |      Chromosome    Start      End  Strand      gene_id          gene_name
          int64  |        category    int64    int64  category    object           object
        -------  ---  ------------  -------  -------  ----------  ---------------  -----------
              0  |               1    11868    14409  +           ENSG00000223972  DDX11L1
              1  |               1    11868    14409  +           ENSG00000223972  DDX11L1
              2  |               1    11868    12227  +           ENSG00000223972  DDX11L1
        PyRanges with 3 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        When requesting non-existing chromosome or strand or ranges an empty PyRanges is returned:

        >>> gr.loci["3"]
        index    |    Chromosome    Start    End      Strand      gene_id    gene_name
        int64    |    category      int64    int64    category    object     object
        -------  ---  ------------  -------  -------  ----------  ---------  -----------
        PyRanges with 0 rows, 6 columns, and 1 index columns.
        Contains 0 chromosomes and 0 strands.

        >>> gr2 = pr.PyRanges({"Chromosome": ["chr1", "chr2"], "Start": [1, 2], "End": [4, 5],
        ...                    "Strand": [".", "+"], "Score":[10, 12], "Id":["a", "b"]})
        >>> gr2.loci["chr2"]
          index  |    Chromosome      Start      End  Strand      Score  Id
          int64  |    object          int64    int64  object      int64  object
        -------  ---  ------------  -------  -------  --------  -------  --------
              1  |    chr2                2        5  +              12  b
        PyRanges with 1 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        The loci operator can also be used for assignment, using a same-sized PyRanges:

        >>> gr2.loci["chr2"] = gr2.loci["chr2"].copy().assign(Chromosome="xxx")
        >>> gr2
          index  |    Chromosome      Start      End  Strand      Score  Id
          int64  |    object          int64    int64  object      int64  object
        -------  ---  ------------  -------  -------  --------  -------  --------
              0  |    chr1                1        4  .              10  a
              1  |    xxx                 2        5  +              12  b
        PyRanges with 2 rows, 6 columns, and 1 index columns.
        Contains 2 chromosomes and 2 strands (including non-genomic strands: .).

        For more flexible assignment, you can employ Pandas loc using the index of the loci output:

        >>> c = gr2.loci["chr1"]
        >>> gr2.loc[c.index, "Score"] = 100
        >>> gr2
          index  |    Chromosome      Start      End  Strand      Score  Id
          int64  |    object          int64    int64  object      int64  object
        -------  ---  ------------  -------  -------  --------  -------  --------
              0  |    chr1                1        4  .             100  a
              1  |    xxx                 2        5  +              12  b
        PyRanges with 2 rows, 6 columns, and 1 index columns.
        Contains 2 chromosomes and 2 strands (including non-genomic strands: .).

        When providing only strand, or strand and a slice, only valid genomic strands (i.e. '+', '-') are searched for:

        >>> gr2.loci['+']
          index  |    Chromosome      Start      End  Strand      Score  Id
          int64  |    object          int64    int64  object      int64  object
        -------  ---  ------------  -------  -------  --------  -------  --------
              1  |    xxx                 2        5  +              12  b
        PyRanges with 1 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        >>> gr2.loci['.']
        index    |    Chromosome    Start    End      Strand    Score    Id
        int64    |    object        int64    int64    object    int64    object
        -------  ---  ------------  -------  -------  --------  -------  --------
        PyRanges with 0 rows, 6 columns, and 1 index columns.
        Contains 0 chromosomes and 0 strands.

        You can use None to match all values, useful to force the non-ambiguous syntax that can match non-genomic strands:

        >>> gr2.loci[None, '.']
          index  |    Chromosome      Start      End  Strand      Score  Id
          int64  |    object          int64    int64  object      int64  object
        -------  ---  ------------  -------  -------  --------  -------  --------
              0  |    chr1                1        4  .             100  a
        PyRanges with 1 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands (including non-genomic strands: .).

        Do not try to use loci to access columns: the key is interpreted as a chromosome, resulting in empty output:

        >>> gr2.loci["Score"]
        index    |    Chromosome    Start    End      Strand    Score    Id
        int64    |    object        int64    int64    object    int64    object
        -------  ---  ------------  -------  -------  --------  -------  --------
        PyRanges with 0 rows, 6 columns, and 1 index columns.
        Contains 0 chromosomes and 0 strands.

        >>> gr2.loci[["Score", "Id"]]
        Traceback (most recent call last):
        ...
        TypeError: The loci accessor does not accept a list. If you meant to retrieve columns, use get_with_loc_columns instead.

        """
        return self._loci

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

        >>> pr.options.get_option('max_rows_to_show')
        8
        >>> pr.options.get_option('console_width')
        120

        >>> gr2 = gr.copy()
        >>> gr2.loc[:, "Strand"] = ["+", "-", "X"]
        >>> gr = pr.concat([gr2, gr, gr])
        >>> gr  # If a PyRanges has more than eight rows the repr is truncated in the middle.
        index    |    Chromosome    Start    End      Strand
        int64    |    int64         int64    int64    object
        -------  ---  ------------  -------  -------  --------
        0        |    3             0        10       +
        1        |    2             100      125      -
        2        |    1             250      251      X
        0        |    3             0        10       .
        ...      |    ...           ...      ...      ...
        2        |    1             250      251      /
        0        |    3             0        10       .
        1        |    2             100      125      ^
        2        |    1             250      251      /
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

    def boundaries(
        self,
        transcript_id: VALID_BY_TYPES = None,
        use_strand: VALID_USE_STRAND_TYPE = "auto",
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

        use_strand: {"auto", True, False}, default: "auto"
            Whether to cluster only intervals on the same strand.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        Returns
        -------
        PyRanges
            One interval per group, with the min(Start) and max(End) of the group

        See Also
        --------
        PyRanges.complement : return the internal complement of intervals, i.e. its introns.

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

        """
        from pyranges.methods.boundaries import _bounds

        by = prepare_by_single(
            self,
            use_strand=use_strand,
            match_by=transcript_id,
        )

        by = [*dict.fromkeys(by).keys()]

        return mypy_ensure_pyranges(
            _bounds(
                df=self,
                by=by,
            )
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
            Length by which the criteria of overlap are loosened.
            A value of 1 clusters also bookended intervals.
            Higher slack values cluster more distant intervals (with a maximum distance of slack-1 between them).

        cluster_column:
            Name the cluster column added in output. Default: "Cluster"

        Returns
        -------
        PyRanges
            PyRanges with an ID-column "Cluster" added.

        See Also
        --------
        PyRanges.merge: combine overlapping intervals into one

        Examples
        --------
        >>> gr = pr.PyRanges(dict(Chromosome=1, Start=[5, 6, 12, 16, 20, 22, 24], End=[9, 8, 16, 18, 23, 25, 27]))
        >>> gr
          index  |      Chromosome    Start      End
          int64  |           int64    int64    int64
        -------  ---  ------------  -------  -------
              0  |               1        5        9
              1  |               1        6        8
              2  |               1       12       16
              3  |               1       16       18
              4  |               1       20       23
              5  |               1       22       25
              6  |               1       24       27
        PyRanges with 7 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.cluster()
          index  |      Chromosome    Start      End    Cluster
          int64  |           int64    int64    int64     uint32
        -------  ---  ------------  -------  -------  ---------
              0  |               1        5        9          0
              1  |               1        6        8          0
              2  |               1       12       16          1
              3  |               1       16       18          2
              4  |               1       20       23          3
              5  |               1       22       25          3
              6  |               1       24       27          3
        PyRanges with 7 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        Slack=1 will cluster also bookended intervals:

        >>> gr.cluster(slack=1)
          index  |      Chromosome    Start      End    Cluster
          int64  |           int64    int64    int64     uint32
        -------  ---  ------------  -------  -------  ---------
              0  |               1        5        9          0
              1  |               1        6        8          0
              2  |               1       12       16          1
              3  |               1       16       18          1
              4  |               1       20       23          2
              5  |               1       22       25          2
              6  |               1       24       27          2
        PyRanges with 7 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        Higher values of slack will cluster more distant intervals:

        >>> gr.cluster(slack=3)
          index  |      Chromosome    Start      End    Cluster
          int64  |           int64    int64    int64     uint32
        -------  ---  ------------  -------  -------  ---------
              0  |               1        5        9          0
              1  |               1        6        8          0
              2  |               1       12       16          1
              3  |               1       16       18          1
              4  |               1       20       23          1
              5  |               1       22       25          1
              6  |               1       24       27          1
        PyRanges with 7 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        result = super().cluster(
            cluster_column=cluster_column,
            match_by=prepare_by_single(self, use_strand=use_strand, match_by=match_by),
            slack=slack,
        )

        return mypy_ensure_pyranges(result)

    def copy(self, *args, **kwargs) -> "pr.PyRanges":
        """Return a copy of the PyRanges."""
        return mypy_ensure_pyranges(super().copy(*args, **kwargs))

    def _count_overlaps(
        self,
        other: "PyRanges",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        *,
        match_by: str | list[str] | None = None,
        slack: int = 0,
    ) -> pd.Series:
        strand_behavior = validate_and_convert_strand_behavior(self, other, strand_behavior)
        _other, by = prepare_by_binary(self, other=other, strand_behavior=strand_behavior, match_by=match_by)
        return super().count_overlaps(_other, match_by=by, slack=slack)

    def count_overlaps(  # type: ignore[override]
        self,
        other: "PyRanges",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        *,
        match_by: str | list[str] | None = None,
        slack: int = 0,
        overlap_col: str = "Count",
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

        slack : int, default 0
            Temporarily lengthen intervals in self before searching for overlaps.

        keep_nonoverlapping : bool, default True
            Keep intervals without overlaps.

        overlap_col : str, default "Count"
            Name of column with overlap counts.

        Returns
        -------
        PyRanges
            PyRanges with a column of overlaps added.

        See Also
        --------
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

        >>> f1.count_overlaps(f2)
          index  |    Chromosome      Start      End  Strand         Count
          int64  |    category        int64    int64  category      uint32
        -------  ---  ------------  -------  -------  ----------  --------
              0  |    chr1                3        6  +                  0
              1  |    chr1                5        7  -                  1
              2  |    chr1                8        9  +                  0
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> f1.count_overlaps(f2, slack=1, strand_behavior="ignore")
          index  |    Chromosome      Start      End  Strand         Count
          int64  |    category        int64    int64  category      uint32
        -------  ---  ------------  -------  -------  ----------  --------
              0  |    chr1                3        6  +                  1
              1  |    chr1                5        7  -                  1
              2  |    chr1                8        9  +                  0
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> annotation = pr.example_data.ensembl_gtf.get_with_loc_columns(['transcript_id', 'Feature'])
        >>> reads = pr.random(1000, chromsizes={'1':150000}, strand=False, seed=123)
        >>> annotation.count_overlaps(reads, overlap_col="NumberOverlaps")
        index    |    Chromosome    Start    End      Strand      transcript_id    Feature     NumberOverlaps
        int64    |    category      int64    int64    category    object           category    uint32
        -------  ---  ------------  -------  -------  ----------  ---------------  ----------  ----------------
        0        |    1             11868    14409    +           nan              gene        17
        1        |    1             11868    14409    +           ENST00000456328  transcript  17
        2        |    1             11868    12227    +           ENST00000456328  exon        3
        3        |    1             12612    12721    +           ENST00000456328  exon        1
        ...      |    ...           ...      ...      ...         ...              ...         ...
        7        |    1             120724   133723   -           ENST00000610542  transcript  76
        8        |    1             133373   133723   -           ENST00000610542  exon        1
        9        |    1             129054   129223   -           ENST00000610542  exon        3
        10       |    1             120873   120932   -           ENST00000610542  exon        1
        PyRanges with 11 rows, 7 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        """
        _self = self.copy()
        _self.loc[:, overlap_col] = _self._count_overlaps(  # noqa: SLF001
            other, strand_behavior=strand_behavior, match_by=match_by, slack=slack
        )
        return _self

    # to do: optimize, doesn't need to split by chromosome, only strand and only if ext_3/5
    def extend(
        self,
        ext: int | None = None,
        ext_3: int | None = None,
        ext_5: int | None = None,
        transcript_id: VALID_BY_TYPES = None,
        use_strand: VALID_USE_STRAND_TYPE = "auto",
    ) -> "PyRanges":
        """Extend the intervals from the 5' and/or 3' ends.

        The Strand (if valid) is considered when extending the intervals:
        a 5' extension applies to the Start of a "+" strand interval and to the End of a "-" strand interval.

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

        use_strand: {"auto", True, False}, default: "auto"
            If False, ignore strand information when extending intervals.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        See Also
        --------
        PyRanges.subsequence : obtain subsequences of intervals
        PyRanges.spliced_subsequence : obtain subsequences of intervals, providing transcript-level coordinates
        PyRanges.five_end : return the 5' end of intervals or transcripts
        PyRanges.three_end : return the 3' end of intervals or transcripts

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
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                4        5  +
              1  |    chr1                9        8  +
              2  |    chr1                6        6  -
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.
        Invalid ranges:
          * 2 intervals are empty or negative length (end <= start). See indexes: 1, 2

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
        import ruranges

        if ext is not None == (ext_3 is not None or ext_5 is not None):
            msg = "Must use at least one and not both of ext and ext3 or ext5."
            raise ValueError(msg)

        use_strand = validate_and_convert_use_strand(self, use_strand) if (ext_3 or ext_5) else False

        groups = factorize(self, transcript_id) if transcript_id is not None else np.arange(len(self), dtype=np.uint32)

        starts, ends = ruranges.extend(
            groups=groups,  # type: ignore[arg-type]
            starts=self[START_COL].to_numpy(),
            ends=self[END_COL].to_numpy(),
            negative_strand=(self[STRAND_COL] == REVERSE_STRAND).to_numpy(),
            ext=ext,
            ext_3=ext_3,
            ext_5=ext_5,
        )

        result = self.copy()
        result.loc[:, START_COL] = starts
        result.loc[:, END_COL] = ends
        return mypy_ensure_pyranges(result)

    def five_end(
        self,
        transcript_id: VALID_BY_TYPES = None,
        slack: int = 0,
    ) -> "PyRanges":
        """Return the five prime end of intervals.

        The five prime end is the start of a forward strand or the end of a reverse strand.

        Parameters
        ----------
        transcript_id : str or list of str, default: None
            Optional column name(s). If provided, the five prime end is calculated for each
            group of intervals.

        slack : int, default 0
            Lengthen (or contract) intervals.

        Returns
        -------
        PyRanges

            PyRanges with the five prime ends

        Note
        ----
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
              2  |    chr1                6        7  -         b
        PyRanges with 2 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        """
        if not self.strand_valid:
            msg = f"Need PyRanges with valid strands ({VALID_GENOMIC_STRAND_INFO}) to find 5'."
            raise AssertionError(msg)

        return mypy_ensure_pyranges(
            process_site(self, slack=slack, site="tss")
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

    def join_ranges(  # type: ignore[override]
        self,
        other: "PyRanges",
        *,
        multiple: VALID_OVERLAP_TYPE = "all",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        join_type: VALID_JOIN_TYPE = "inner",
        match_by: VALID_BY_TYPES = None,
        contained_intervals_only: bool = False,
        slack: int = 0,
        suffix: str = JOIN_SUFFIX,
        report_overlap_column: str | None = None,
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

        multiple : {"all", "first", "last"}, default "all"
            What intervals to report when multiple intervals in 'other' overlap with the same interval in self.
            The default "all" reports all overlapping subintervals, which will have duplicate indices.
            "first" reports only, for each interval in self, the overlapping subinterval with smallest Start in 'other'
            "last" reports only the overlapping subinterval with the biggest End in 'other'

        contained_intervals_only : bool, default False
            Whether to report only intervals that are entirely contained in an interval of 'other'.

        match_by : str or list, default None
            If provided, only intervals with an equal value in column(s) `match_by` may be joined.

        report_overlap_column : str or None
            Report amount of overlap in base pairs using column name

        slack : int, default 0
            Before joining, temporarily extend intervals in self by this much on both ends.

        suffix : str or tuple, default "_b"
            Suffix to give overlapping columns in other.

        Returns
        -------
        PyRanges

            A PyRanges appended with columns of another.

        Note
        ----
        The indices of the two input PyRanges are not preserved in output.
        The chromosome column from other will never be reported as it is always the same as in self.
        Whether the strand column from other is reported depends on the strand_behavior.


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

        Note that since some start and end columns are nan, a regular DataFrame is returned.

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
              0  |    chr1                1        3  a         chr1                    2        3  a
              0  |    chr1                1        3  a         chr1                    2        9  b
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

        With slack 1, bookended features are joined (see row 1):

        >>> f1.join_ranges(f2, slack=1)
          index  |    Chromosome      Start      End  Name       Chromosome_b      Start_b    End_b  Name_b
          int64  |    object          int64    int64  object     object              int64    int64  object
        -------  ---  ------------  -------  -------  ---------  --------------  ---------  -------  --------
              0  |    chr1                3        6  interval1  chr1                    6        7  b
              2  |    chr1                5        7  interval2  chr1                    6        7  b
        PyRanges with 2 rows, 8 columns, and 1 index columns.
        Contains 1 chromosomes.


        >>> f1.join_ranges(f2, report_overlap_column="Overlap")
          index  |    Chromosome      Start      End  Name       Chromosome_b      Start_b    End_b  Name_b      Overlap
          int64  |    object          int64    int64  object     object              int64    int64  object        int64
        -------  ---  ------------  -------  -------  ---------  --------------  ---------  -------  --------  ---------
              2  |    chr1                5        7  interval2  chr1                    6        7  b                 1
        PyRanges with 1 rows, 9 columns, and 1 index columns.
        Contains 1 chromosomes.

        Allowing slack in overlaps may result in 0 or negative Overlap values:

        >>> f1.join_ranges(f2, report_overlap_column="Overlap", slack=2)
          index  |    Chromosome      Start      End  Name       Chromosome_b      Start_b    End_b  Name_b      Overlap
          int64  |    object          int64    int64  object     object              int64    int64  object        int64
        -------  ---  ------------  -------  -------  ---------  --------------  ---------  -------  --------  ---------
              0  |    chr1                3        6  interval1  chr1                    1        2  a                -1
              0  |    chr1                3        6  interval1  chr1                    6        7  b                 0
              1  |    chr1                8        9  interval3  chr1                    6        7  b                -1
              2  |    chr1                5        7  interval2  chr1                    6        7  b                 1
        PyRanges with 4 rows, 9 columns, and 1 index columns (with 1 index duplicates).
        Contains 1 chromosomes.

        """
        _other, by = prepare_by_binary(self, other=other, strand_behavior=strand_behavior, match_by=match_by)
        gr = super().join_ranges(
            other=_other,
            match_by=by,
            join_type=join_type,
            multiple=multiple,
            contained_intervals_only=contained_intervals_only,
            slack=slack,
            suffix=suffix,
            report_overlap_column=report_overlap_column,
        )

        return mypy_ensure_pyranges(gr)

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
        PyRanges.length : return the total length of all intervals combined

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

    def map_to_global(
        self,
        gr: "PyRanges",
        global_on: str,
        *,
        local_on: str = "Chromosome",
        use_strand: VALID_USE_STRAND_TYPE = USE_STRAND_DEFAULT,
    ) -> "PyRanges":
        """Map intervals from a *local* reference frame (e.g. transcript) to *global* coordinates (e.g. genomic).

        The self PyRanges object is *local* in the sense that its
        ``Chromosome`` column stores an **identifier** (e.g. a
        transcript ID), so that its interval coordinates are expressed
        relative to that identifier.
        The *global* object *gr* supplies the absolute genomic coordinates
        of every interval group (e.g. transcripts, potentially with multiple
        exons). The function returns the self intervals in the reference
        system of the global object.

        The strand of returned intervals is the product of the strand of
        the corresponding local and global intervals (e.g. +/- => -)

        Unused rows in *gr* (identifiers never referenced by ``self``) are ignored.

        Parameters
        ----------
        gr : PyRanges
            Intervals in global reference system (e.g. transcript annotation
            in genomic coordinates).
        global_on : str
            Column in *gr* that holds the identifiers contained in
            ``self.Chromosome``.
        local_on : str, default "Chromosome"
            Column in `self` that holds the identifier to be lifted. Change this
            if your identifiers live in a different column.
        use_strand: {"auto", True, False}, default: "auto"
            Only map when strands in local and global match.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        Returns
        -------
        PyRanges
            Intervals in genomic coordinates, maintaining order, index, and
            metadata columns of self.

        Warning
        -------
        A single local interval will give rise to multiple intervals in output
        if it overlaps discontinuities (i.e. introns) in global coordinates.
        This will generate duplicated indices. To avoid them,
        run pandas dataframe method ``reset_index`` on the output.

        Examples
        --------
        >>> gr = pr.PyRanges(pd.DataFrame({
        ...     "Chromosome": ["chr1","chr1","chr1","chr1"],
        ...     "Start": [100, 300, 1000, 1100],
        ...     "End": [200, 400, 1050, 1200],
        ...     "Strand": ["+","+", "-", "-"],
        ...     "transcript_id": ["tx1","tx1","tx2","tx2"],
        ... }))
        >>> tr = pr.PyRanges(pd.DataFrame({
        ...     "Chromosome": ["tx1","tx1","tx1","tx2","tx2"],
        ...     "Start": [0, 120, 160, 0, 100],
        ...     "End": [80, 140, 170, 20, 130],
        ...     "Strand": ["-","-", "+", "+", "+"],
        ...     "label": ["a","b","c","d","e"],
        ... }))
        >>> gr
          index  |    Chromosome      Start      End  Strand    transcript_id
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  ---------------
              0  |    chr1              100      200  +         tx1
              1  |    chr1              300      400  +         tx1
              2  |    chr1             1000     1050  -         tx2
              3  |    chr1             1100     1200  -         tx2
        PyRanges with 4 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> tr
          index  |    Chromosome      Start      End  Strand    label
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  --------
              0  |    tx1                 0       80  -         a
              1  |    tx1               120      140  -         b
              2  |    tx1               160      170  +         c
              3  |    tx2                 0       20  +         d
              4  |    tx2               100      130  +         e
        PyRanges with 5 rows, 5 columns, and 1 index columns.
        Contains 2 chromosomes and 2 strands.

        >>> tr.map_to_global(gr, "transcript_id")
          index  |    Chromosome      Start      End  Strand    label
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  --------
              0  |    chr1              100      180  -         a
              1  |    chr1              320      340  -         b
              2  |    chr1              360      370  +         c
              3  |    chr1             1180     1200  -         d
              4  |    chr1             1020     1050  -         e
        PyRanges with 5 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        Extra columns are preserved:

        >>> tr.assign(tag=7).map_to_global(gr, "transcript_id").tag.unique()
        array([7])

        A local interval that spans an exon junction is split; its index is
        duplicated in the output.

        >>> tr2 = pr.PyRanges(pd.DataFrame({
        ...     "Chromosome":["tx1","tx2","tx2"],
        ...     "Start": [90, 80, 50],
        ...     "End": [110, 120, 120],
        ...     "Strand": ["+","+", "-"],
        ...     "label": ["q","w","e"],
        ... }))
        >>> tr2.map_to_global(gr, "transcript_id")    # doctest: +ELLIPSIS
          index  |    Chromosome      Start      End  Strand    label
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  --------
              0  |    chr1              190      200  +         q
              0  |    chr1              300      310  +         q
              1  |    chr1             1100     1120  -         w
              1  |    chr1             1030     1050  -         w
              2  |    chr1             1100     1150  +         e
              2  |    chr1             1030     1050  +         e
        PyRanges with 6 rows, 5 columns, and 1 index columns (with 3 index duplicates).
        Contains 1 chromosomes and 2 strands.

        A local interval longer than its transcript is truncated to the
        portion that fits.

        >>> tr3 = pr.PyRanges(pd.DataFrame({
        ...     "Chromosome":["tx1"], "Start":[20], "End":[1000], "Strand":["+"]
        ... }))
        >>> tr3.map_to_global(gr, "transcript_id")
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1              120      200  +
              0  |    chr1              300      400  +
        PyRanges with 2 rows, 4 columns, and 1 index columns (with 1 index duplicates).
        Contains 1 chromosomes and 1 strands.

        """
        return map_to_global(
            local_gr=self,
            global_gr=gr,
            global_on=global_on,
            local_on=local_on,
            use_strand=use_strand and STRAND_COL in self,
        )

    def max_disjoint(
        self,
        use_strand: VALID_USE_STRAND_TYPE = "auto",
        *,
        slack: int = 0,
        match_by: VALID_BY_TYPES = None,
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
            Length by which the criteria of overlap are loosened.
            A value of 1 implies that bookended intervals are considered overlapping.
            Higher slack values allow more distant intervals (with a maximum distance of slack-1 between them).

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
        PyRanges.cluster : annotate overlapping intervals with common ID

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
        use_strand = validate_and_convert_use_strand(self, use_strand)

        result = super().max_disjoint(
            match_by=prepare_by_single(self, use_strand=use_strand, match_by=match_by),
            slack=slack,
        )
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

        Returns a PyRanges with superintervals that are the union of overlapping intervals.

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
            PyRanges with superintervals. Metadata columns, index, and order are not preserved.

        Note
        ----
        To avoid losing metadata, use cluster instead. If you want to perform a reduction function
        on the metadata, use pandas groupby.

        See Also
        --------
        PyRanges.cluster : annotate overlapping intervals with common ID
        PyRanges.max_disjoint : find the maximal disjoint set of intervals
        PyRanges.split : split intervals into non-overlapping subintervals

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
          index  |      Chromosome    Start      End  Strand         Count
          int64  |        category    int64    int64  category      uint32
        -------  ---  ------------  -------  -------  ----------  --------
              0  |               1    11868    14409  +                  5
              1  |               1   110952   111357  -                  1
              2  |               1   112699   112804  -                  1
              3  |               1   120724   133723  -                  4
        PyRanges with 4 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.merge_overlaps(count_col="Count", match_by="gene_name")
          index  |      Chromosome    Start      End  Strand      gene_name       Count
          int64  |        category    int64    int64  category    object         uint32
        -------  ---  ------------  -------  -------  ----------  -----------  --------
              0  |               1    11868    14409  +           DDX11L1             5
              1  |               1   110952   111357  -           AL627309.1          1
              2  |               1   112699   112804  -           AL627309.1          1
              3  |               1   120724   133723  -           AL627309.1          4
        PyRanges with 4 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.


        """
        result = super().merge_overlaps(
            count_col=count_col,
            match_by=prepare_by_single(self, use_strand=use_strand, match_by=match_by),
            slack=slack,
        )

        return mypy_ensure_pyranges(result)

    def nearest(  # type: ignore[override]
        self,
        other: "PyRanges",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        direction: VALID_NEAREST_TYPE = "any",
        *,
        k: int = 1,
        match_by: VALID_BY_TYPES = None,
        suffix: str = JOIN_SUFFIX,
        exclude_overlaps: bool = False,
        dist_col: str | None = "Distance",
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

        exclude_overlaps : bool, default True
            Whether to not report intervals of others that overlap with self as the nearest ones.

        direction : {"any", "upstream", "downstream"}, default "any", i.e. both directions
            Whether to only look for nearest in one direction.

        match_by : str or list, default None
            If provided, only intervals with an equal value in column(s) `match_by` may be matched.

        k : int, default 1
            Number of nearest intervals to fetch.

        suffix : str, default "_b"
            Suffix to give columns with shared name in other.

        dist_col : str or None
            Optional column to store the distance in.

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

        >>> f2 = pr.PyRanges(dict(Chromosome="chr1", Start=[1, 6, 20], End=[2, 7, 22], Strand=["+", "-", "+"]))
        >>> f2
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        2  +
              1  |    chr1                6        7  -
              2  |    chr1               20       22  +
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> f1.nearest(f2)
          index  |    Chromosome      Start      End  Strand      Chromosome_b      Start_b    End_b  Strand_b      Distance
          int64  |    category        int64    int64  category    object              int64    int64  object           int64
        -------  ---  ------------  -------  -------  ----------  --------------  ---------  -------  ----------  ----------
              0  |    chr1                3        6  +           chr1                    1        2  +                    2
              1  |    chr1                5        7  -           chr1                    6        7  -                    0
              2  |    chr1                8        9  +           chr1                    1        2  +                    7
        PyRanges with 3 rows, 9 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> f1.nearest(f2, strand_behavior='ignore')
          index  |    Chromosome      Start      End  Strand      Chromosome_b      Start_b    End_b  Strand_b      Distance
          int64  |    category        int64    int64  category    object              int64    int64  object           int64
        -------  ---  ------------  -------  -------  ----------  --------------  ---------  -------  ----------  ----------
              0  |    chr1                3        6  +           chr1                    6        7  -                    1
              1  |    chr1                5        7  -           chr1                    6        7  -                    0
              2  |    chr1                8        9  +           chr1                    6        7  -                    2
        PyRanges with 3 rows, 9 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> f1.nearest(f2, k=2, strand_behavior='ignore')
          index  |    Chromosome      Start      End  Strand      Chromosome_b      Start_b    End_b  Strand_b      Distance
          int64  |    category        int64    int64  category    object              int64    int64  object           int64
        -------  ---  ------------  -------  -------  ----------  --------------  ---------  -------  ----------  ----------
              0  |    chr1                3        6  +           chr1                    6        7  -                    1
              0  |    chr1                3        6  +           chr1                    1        2  +                    2
              1  |    chr1                5        7  -           chr1                    6        7  -                    0
              1  |    chr1                5        7  -           chr1                    1        2  +                    4
              2  |    chr1                8        9  +           chr1                    6        7  -                    2
              2  |    chr1                8        9  +           chr1                    1        2  +                    7
        PyRanges with 6 rows, 9 columns, and 1 index columns (with 3 index duplicates).
        Contains 1 chromosomes and 2 strands.

        >>> f1.nearest(f2, strand_behavior='ignore', exclude_overlaps=True)
          index  |    Chromosome      Start      End  Strand      Chromosome_b      Start_b    End_b  Strand_b      Distance
          int64  |    category        int64    int64  category    object              int64    int64  object           int64
        -------  ---  ------------  -------  -------  ----------  --------------  ---------  -------  ----------  ----------
              0  |    chr1                3        6  +           chr1                    6        7  -                    1
              1  |    chr1                5        7  -           chr1                    1        2  +                    4
              2  |    chr1                8        9  +           chr1                    6        7  -                    2
        PyRanges with 3 rows, 9 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> f1.nearest(f2, direction='downstream')
          index  |    Chromosome      Start      End  Strand      Chromosome_b      Start_b    End_b  Strand_b      Distance
          int64  |    category        int64    int64  category    object              int64    int64  object           int64
        -------  ---  ------------  -------  -------  ----------  --------------  ---------  -------  ----------  ----------
              0  |    chr1                3        6  +           chr1                   20       22  +                   15
              2  |    chr1                8        9  +           chr1                   20       22  +                   12
              1  |    chr1                5        7  -           chr1                    6        7  -                    0
        PyRanges with 3 rows, 9 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        If an input interval has no suitable nearest interval, these rows are dropped:

        >>> f1.nearest(f2, direction='upstream', exclude_overlaps=True)
          index  |    Chromosome      Start      End  Strand      Chromosome_b      Start_b    End_b  Strand_b      Distance
          int64  |    category        int64    int64  category    object              int64    int64  object           int64
        -------  ---  ------------  -------  -------  ----------  --------------  ---------  -------  ----------  ----------
              0  |    chr1                3        6  +           chr1                    1        2  +                    2
              2  |    chr1                8        9  +           chr1                    1        2  +                    7
        PyRanges with 2 rows, 9 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        """
        _other, by = prepare_by_binary(self, other=other, strand_behavior=strand_behavior, match_by=match_by)

        if direction == NEAREST_ANY_DIRECTION:
            res = super().nearest(
                other=_other,
                match_by=by,
                suffix=suffix,
                exclude_overlaps=exclude_overlaps,
                k=k,
                dist_col=dist_col,
                direction="any",
            )
            return mypy_ensure_pyranges(res)

        fwd_self, rev_self = split_on_strand(self)
        if direction == NEAREST_DOWNSTREAM:
            res = RangeFrame(fwd_self).nearest(
                other=_other,
                match_by=by,
                suffix=suffix,
                exclude_overlaps=exclude_overlaps,
                k=k,
                dist_col=dist_col,
                direction="forward",
            )
            res2 = RangeFrame(rev_self).nearest(
                other=_other,
                match_by=by,
                suffix=suffix,
                exclude_overlaps=exclude_overlaps,
                k=k,
                dist_col=dist_col,
                direction="forward",
            )
        elif direction == NEAREST_UPSTREAM:
            res = RangeFrame(fwd_self).nearest(
                other=RangeFrame(_other),
                match_by=by,
                suffix=suffix,
                exclude_overlaps=exclude_overlaps,
                k=k,
                dist_col=dist_col,
                direction="backward",
            )
            res2 = RangeFrame(rev_self).nearest(
                other=RangeFrame(_other),
                match_by=by,
                suffix=suffix,
                exclude_overlaps=exclude_overlaps,
                k=k,
                dist_col=dist_col,
                direction="backward",
            )
        else:
            msg = f"Invalid direction: {direction}"
            raise ValueError(msg)

        return mypy_ensure_pyranges(
            pd.concat(
                [res, res2],
            )
        )

    def overlap(  # type: ignore[override]
        self,
        other: "PyRanges",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        multiple: VALID_OVERLAP_TYPE = "all",
        slack: int = 0,
        *,
        contained_intervals_only: bool = False,
        match_by: VALID_BY_TYPES = None,
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

        multiple : {"all", "first", "last"}, default "all"
            What intervals to report when multiple intervals in 'other' overlap with the same interval in self.
            The default "all" reports all overlapping subintervals, which will have duplicate indices.
            "first" reports only, for each interval in self, the overlapping subinterval with smallest Start in 'other'
            "last" reports only the overlapping subinterval with the biggest End in 'other'

        slack : int, default 0
            Intervals in self are temporarily extended by slack on both ends before overlap is calculated, so that
            we allow non-overlapping intervals to be considered overlapping if they are within less than slack distance
            e.g. slack=1 reports bookended intervals.

        contained_intervals_only : bool, default False
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
        >>> gr
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        3  A
              1  |    chr1                1        3  a
              2  |    chr2                4        9  b
              3  |    chr1               10       11  c
              4  |    chr3                0        1  d
        PyRanges with 5 rows, 4 columns, and 1 index columns.
        Contains 3 chromosomes.

        >>> gr2
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                2        3
              1  |    chr1                2        9
              2  |    chr2                1       10
        PyRanges with 3 rows, 3 columns, and 1 index columns.
        Contains 2 chromosomes.

        >>> gr.overlap(gr2)
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        3  A
              0  |    chr1                1        3  A
              1  |    chr1                1        3  a
              1  |    chr1                1        3  a
              2  |    chr2                4        9  b
        PyRanges with 5 rows, 4 columns, and 1 index columns (with 2 index duplicates).
        Contains 2 chromosomes.

        >>> gr.overlap(gr2, slack=2, multiple="first")
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        3  A
              1  |    chr1                1        3  a
              2  |    chr2                4        9  b
              3  |    chr1               10       11  c
        PyRanges with 4 rows, 4 columns, and 1 index columns.
        Contains 2 chromosomes.

        >>> gr.overlap(gr2, multiple="last")
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        3  A
              1  |    chr1                1        3  a
              2  |    chr2                4        9  b
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 2 chromosomes.

        >>> gr.overlap(gr2, contained_intervals_only=True)
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              2  |    chr2                4        9  b
        PyRanges with 1 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.overlap(gr2, contained_intervals_only=True, slack=-2)
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        3  A
              1  |    chr1                1        3  a
              2  |    chr2                4        9  b
        PyRanges with 3 rows, 4 columns, and 1 index columns.
        Contains 2 chromosomes.

        >>> gr.overlap(gr2, contained_intervals_only=True, multiple="first")
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              2  |    chr2                4        9  b
        PyRanges with 1 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

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
        _other, by = prepare_by_binary(self, other=other, strand_behavior=strand_behavior, match_by=match_by)

        gr = super().overlap(
            _other,
            match_by=by,
            slack=slack,
            multiple=multiple,
            contained_intervals_only=contained_intervals_only,
        )

        return mypy_ensure_pyranges(gr)

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
            The default "all" reports all overlapping subintervals.
            "first" reports only, for each merged self interval, the overlapping 'other' subinterval with smallest Start
            "last" reports only the overlapping subinterval with the biggest End in 'other'

        Returns
        -------
        PyRanges
            A PyRanges with overlapping subintervals. Input index is not preserved.
            No columns other than Chromosome, Start, End, and optionally Strand are returned.

        See Also
        --------
        PyRanges.set_union : set-theoretical union
        PyRanges.intersect : find overlapping subintervals
        PyRanges.overlap : report overlapping intervals
        PyRanges.merge_overlaps : merge overlapping intervals

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

        >>> r1.set_intersect(r2, multiple='first')
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                7        9
              1  |    chr1               20       22
        PyRanges with 2 rows, 3 columns, and 1 index columns.
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

        """
        strand_behavior = validate_and_convert_strand_behavior(self, other, strand_behavior)

        use_strand = use_strand_from_validated_strand_behavior(self, other, strand_behavior)
        self_clusters = self.merge_overlaps(use_strand=use_strand and self.has_strand)
        other_clusters = other.merge_overlaps(use_strand=use_strand and other.has_strand)
        result = self_clusters.intersect(other_clusters, strand_behavior=strand_behavior, multiple=multiple)
        return mypy_ensure_pyranges(result.reset_index(drop=True))

    def set_union(self, other: "PyRanges", strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto") -> "PyRanges":
        """Return set-theoretical union.

        Returns the regions present in either self or other.
        Both PyRanges are merged first.

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
            A PyRanges with the union of intervals. Input index is not preserved.
            No columns other than Chromosome, Start, End, and optionally Strand are returned.

        See Also
        --------
        PyRanges.set_intersect : set-theoretical intersection
        PyRanges.overlap : report overlapping intervals
        PyRanges.merge_overlaps : merge overlapping intervals

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
              0  |    chr1                1        9
              1  |    chr1                9       10
              2  |    chr1               10       11
        PyRanges with 3 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        Merging bookended intervals:

        >>> gr.set_union(gr2).merge_overlaps(slack=1)
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                1       11
        PyRanges with 1 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.


        """
        if self.empty and other.empty:
            return mypy_ensure_pyranges(self.copy())

        strand_behavior = validate_and_convert_strand_behavior(self, other, strand_behavior)
        use_strand = use_strand_from_validated_strand_behavior(self, other, strand_behavior)
        _self = self
        if not use_strand:
            _self = _self.remove_strand()
            other = other.remove_strand()
        elif strand_behavior == STRAND_BEHAVIOR_OPPOSITE:
            other = mypy_ensure_pyranges(
                other.assign(
                    **{
                        STRAND_COL: other[STRAND_COL].replace(
                            {FORWARD_STRAND: REVERSE_STRAND, REVERSE_STRAND: FORWARD_STRAND},
                        ),
                    },
                ),
            )

        gr = pr.concat([_self, other])

        return gr.merge_overlaps(use_strand=use_strand)

    def sort_ranges(  # type: ignore[override]
        self,
        match_by: VALID_BY_TYPES = None,
        *,
        natsort: bool = True,
        use_strand: VALID_USE_STRAND_TYPE = "auto",
    ) -> "PyRanges":
        """Sort PyRanges according to Chromosome, Strand (if present), Start, and End; or by the specified columns.

        If PyRanges is stranded and use_strand is True, intervals on the negative strand are sorted in descending
        order, and End is considered before Start. This is to have a 5' to 3' order.
        For uses not covered by this function, use  DataFrame.sort_values().

        Parameters
        ----------
        match_by : str or list of str, default None
            If provided, sorting occurs by Chromosome, Strand (if present), *by, Start, and End.
            To prioritize columns differently (e.g. Strand before Chromosome), explicitly provide all columns
            in the desired order as part of the 'by' argument.
            You can prepend any column name with '-' to reverse the order of sorting for that column.

        use_strand: {"auto", True, False}, default: "auto"
            Whether negative strand intervals should be sorted in descending order, meaning 5' to 3'.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        natsort : bool, default False
            Whether to use natural sorting for Chromosome column, so that e.g. chr2 < chr11.

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

        >>> p.sort_ranges(natsort=False)
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

        >>> p.sort_ranges(use_strand=False, natsort=False)
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

        >>> p.sort_ranges()
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

        >>> p.sort_ranges(match_by='transcript_id', natsort=False)
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

        >>> res = p.sort_ranges(natsort=False)
        >>> res.sort_values("transcript_id", kind="stable")
          index  |    Chromosome    Strand      Start      End  transcript_id
          int64  |    object        object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              7  |    chr1          +              90      100  t1
              3  |    chr1          -              70       80  t2
              2  |    chr1          -              10       25  t2
              1  |    chr1          +               1       11  t3
              0  |    chr1          +              40       60  t3
              4  |    chr2          +             300      400  t4
              5  |    chr11         +             140      152  t5
              6  |    chr11         +             160      190  t5
        PyRanges with 8 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        Same as before, but 'transcript_id' is sorted in descending order:

        >>> res = p.sort_ranges(natsort=False)
        >>> res.sort_values("transcript_id", kind="stable", ascending=False)
          index  |    Chromosome    Strand      Start      End  transcript_id
          int64  |    object        object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              5  |    chr11         +             140      152  t5
              6  |    chr11         +             160      190  t5
              4  |    chr2          +             300      400  t4
              1  |    chr1          +               1       11  t3
              0  |    chr1          +              40       60  t3
              3  |    chr1          -              70       80  t2
              2  |    chr1          -              10       25  t2
              7  |    chr1          +              90      100  t1
        PyRanges with 8 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        """
        import ruranges

        by = arg_to_list(match_by)

        use_strand = validate_and_convert_use_strand(self, use_strand)

        by = ([CHROM_COL] if STRAND_COL not in self else CHROM_AND_STRAND_COLS) + by

        by_sort_order_as_int = sort_factorize_dict(self, by, use_natsort=natsort)
        idxs = ruranges.sort_intervals(  # type: ignore[attr-defined]
            self[START_COL].to_numpy(),
            self[END_COL].to_numpy(),
            by_sort_order_as_int,
            sort_reverse_direction=None if not use_strand else (self[STRAND_COL] == "-").to_numpy(dtype=bool),
        )
        res = self.take(idxs)  # type: ignore[arg-type]

        return mypy_ensure_pyranges(res)

    def spliced_subsequence(
        self,
        start: int = 0,
        end: int | None = None,
        transcript_id: VALID_BY_TYPES = None,
        use_strand: VALID_USE_STRAND_TYPE = "auto",
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
            If not provided, it is equivalent to PyRanges.subsequence.

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
        PyRanges.window : divide intervals into windows

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

        Get the first 15 nucleotides of *each spliced transcript*, grouping exons by transcript_id:

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

        Get the last 20 nucleotides of each spliced transcript:

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

        Use use_strand=False to treat all intervals as if they were on the + strand:

        >>> p.spliced_subsequence(0, 15, transcript_id='transcript_id', use_strand=False)
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              0  |               1  +               1       11  t1
              1  |               1  +              40       45  t1
              2  |               2  -              10       25  t2
              4  |               3  +             140      152  t3
        PyRanges with 4 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        Get region from 25 to 60 of each spliced transcript, or their existing subportion:

        >>> p.spliced_subsequence(25, 60, transcript_id='transcript_id')
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              1  |               1  +              55       60  t1
        PyRanges with 1 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        Get region of each spliced transcript which excludes their first and last 3 nucleotides:

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
        from pyranges.methods.spliced_subsequence import _spliced_subseq

        by = (
            prepare_by_single(
                self,
                use_strand=use_strand,
                match_by=transcript_id,
            )
            if transcript_id
            else []
        )

        result = _spliced_subseq(
            self,
            by=by,
            force_plus_strand=not use_strand,
            start=start,
            end=end,
        )

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
        PyRanges.merge_overlaps : merge overlapping intervals
        PyRanges.max_disjoint : find the maximal disjoint set of intervals
        PyRanges.split : split intervals into non-overlapping subintervals


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
              4  |    chr1               11       12  -
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
              4  |    chr1               11       12
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
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                3        5  a
              1  |    chr1                5        6  a
              2  |    chr1                6        7  a
              3  |    chr1                5        9  b
              4  |    chr1               11       12  c
        PyRanges with 5 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        import ruranges

        use_strand = validate_and_convert_use_strand(self, use_strand=use_strand)
        by = prepare_by_single(self, use_strand=use_strand, match_by=match_by)
        groups = factorize(self, by)

        idxs, starts, ends = ruranges.split(
            groups=groups,
            starts=self[START_COL].to_numpy(),
            ends=self[END_COL].to_numpy(),
            slack=0,
            between=between,
        )

        res = mypy_ensure_pyranges(self.take(idxs).reset_index(drop=True))  # type: ignore[arg-type]
        if between:
            res = res.remove_nonloc_columns()

        if not use_strand:
            res = res.remove_strand()

        res.loc[:, START_COL] = starts
        res.loc[:, END_COL] = ends

        return mypy_ensure_pyranges(res)

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
        if not self.has_strand:
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
        PyRanges.window : divide intervals into windows

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


        >>> p.subsequence(30, 300, transcript_id='transcript_id')
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              1  |               1  +              40       60  t1
              2  |               2  -               2       13  t2
        PyRanges with 2 rows, 5 columns, and 1 index columns.
        Contains 2 chromosomes and 2 strands.

        >>> p.subsequence(0, 10)
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              0  |               1  +               1       11  t1
              1  |               1  +              40       50  t1
              2  |               2  -               3       13  t2
              3  |               2  -              35       45  t2
              4  |               3  +             140      150  t3
        PyRanges with 5 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        >>> p.subsequence(-5)
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              0  |               1  +              15       20  t1
              1  |               1  +              55       60  t1
              2  |               2  -               2        7  t2
              3  |               2  -              30       35  t2
              4  |               3  +             150      155  t3
        PyRanges with 5 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        Use use_strand=False to treat all intervals as if they were on the + strand:

        >>> p.subsequence(0, 10, use_strand=False)
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

        Get the first 10 nucleotides of *each transcript*, grouping exons by transcript_id:

        >>> p.subsequence(0, 10, transcript_id='transcript_id')
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              0  |               1  +               1       11  t1
              3  |               2  -              35       45  t2
              4  |               3  +             140      150  t3
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        Get the last 20 nucleotides of each transcript:

        >>> p.subsequence(-20, transcript_id='transcript_id')
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              1  |               1  +              40       60  t1
              2  |               2  -               2       13  t2
              4  |               3  +             140      155  t3
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        Get the first 10 nucleotides (at the 5') of *each interval* (each line of the dataframe):

        Get region from 30 to 330 of each transcript, or their existing subportion:

        >>> p.subsequence(30, 300, transcript_id='transcript_id')
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              1  |               1  +              40       60  t1
              2  |               2  -               2       13  t2
        PyRanges with 2 rows, 5 columns, and 1 index columns.
        Contains 2 chromosomes and 2 strands.

        >>> p.subsequence(-20, transcript_id='transcript_id', use_strand=False)
          index  |      Chromosome  Strand      Start      End  transcript_id
          int64  |           int64  object      int64    int64  object
        -------  ---  ------------  --------  -------  -------  ---------------
              1  |               1  +              40       60  t1
              3  |               2  -              30       45  t2
              4  |               3  +             140      155  t3
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 3 chromosomes and 2 strands.

        """
        by = prepare_by_single(self, use_strand=use_strand, match_by=transcript_id) if transcript_id else []

        boundaries = self.boundaries(by, use_strand=use_strand) if by else self
        result = boundaries.spliced_subsequence(
            transcript_id=[],
            use_strand=use_strand,
            start=start,
            end=end,
        )
        return mypy_ensure_pyranges(self.intersect(result))

    def subtract_ranges(  # type: ignore[override]
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

        Returns
        -------
        PyRanges
            PyRanges with subintervals from self that do not overlap with any interval in other.
            Columns and index are preserved.

        Warning
        -------
        The returned Pyranges may have index duplicates. Call .reset_index(drop=True) to fix it.

        See Also
        --------
        PyRanges.overlap : use with invert=True to return all intervals without overlap
        PyRanges.complement : return the internal complement of intervals, i.e. its introns.

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
        _other, by = prepare_by_binary(self, other=other, strand_behavior=strand_behavior, match_by=match_by)

        gr = super().subtract_ranges(
            _other,
            match_by=by,
        )

        return mypy_ensure_pyranges(gr)

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
        use_strand: bool = False,
        overlap_column: str | None = None,
    ) -> "PyRanges":
        """Return overlapping genomic tiles.

        The genome is divided into bookended tiles of length `tile_size`. One tile is returned for each
        interval that overlaps with it, including any metadata from the original intervals.

        Parameters
        ----------
        tile_size : int
            Length of the tiles.

        overlap_column : str, default None
            Name of column to add with the overlap between each bookended tile.

        use_strand: {"auto", True, False}, default: "auto"
            Whether negative strand intervals should be windowed in reverse order.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

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
        pyranges.tile_genome : divide the genome into tiles


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
        int64    |    category      int64    int64    category    category    object       float64
        -------  ---  ------------  -------  -------  ----------  ----------  -----------  -------------
        0        |    1             11800    11900    +           gene        DDX11L1      0.32
        0        |    1             11900    12000    +           gene        DDX11L1      1.0
        0        |    1             12000    12100    +           gene        DDX11L1      1.0
        0        |    1             12100    12200    +           gene        DDX11L1      1.0
        ...      |    ...           ...      ...      ...         ...         ...          ...
        9        |    1             129100   129200   -           exon        AL627309.1   1.0
        9        |    1             129200   129300   -           exon        AL627309.1   0.23
        10       |    1             120800   120900   -           exon        AL627309.1   0.27
        10       |    1             120900   121000   -           exon        AL627309.1   0.32
        PyRanges with 223 rows, 7 columns, and 1 index columns (with 212 index duplicates).
        Contains 1 chromosomes and 2 strands.

        """
        import ruranges

        use_strand = validate_and_convert_use_strand(self, use_strand)

        negative_strand = (self[STRAND_COL] == "-").to_numpy() if use_strand else np.zeros(len(self), dtype=bool)
        indices, starts, ends, overlap_fraction = ruranges.tile(  # type: ignore[attr-defined]
            starts=self[START_COL].to_numpy(),
            ends=self[END_COL].to_numpy(),
            negative_strand=negative_strand,
            tile_size=tile_size,
        )

        res = self.take(indices)  # type: ignore[arg-type]
        res.loc[:, START_COL] = starts
        res.loc[:, END_COL] = ends
        if overlap_column:
            res.loc[:, overlap_column] = overlap_fraction

        # every interval can be processed individually. This may be optimized in the future
        return mypy_ensure_pyranges(res)

    def three_end(
        self,
        transcript_id: VALID_BY_TYPES = None,
        slack: int = 0,
    ) -> "PyRanges":
        """Return the 3'-end.

        The 3'-end is the the end of intervals on the forward strand and the start of intervals on the reverse strand.

        Parameters
        ----------
        transcript_id : str or list of str, default: None
            Optional column name(s). If provided, the three prime end is calculated for each
            group of intervals.

        slack : int, default 0
            Lengthen (or contract) intervals.

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

        >>> gr.three_end(transcript_id='Name')
          index  |    Chromosome      Start      End  Strand    Name
          int64  |    object          int64    int64  object    object
        -------  ---  ------------  -------  -------  --------  --------
              1  |    chr1               13       14  +         a
              2  |    chr1                5        6  -         b
        PyRanges with 2 rows, 5 columns, and 1 index columns.
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

        """
        if not (self.has_strand and self.strand_valid):
            msg = f"Need PyRanges with valid strands ({VALID_GENOMIC_STRAND_INFO}) to find 3'."
            raise AssertionError(msg)

        return mypy_ensure_pyranges(
            process_site(self, slack=slack, site="tes")
            if transcript_id is None
            else self.subsequence(transcript_id=transcript_id, start=-1, use_strand=True)
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

        File contents::

            chr1	1	5	.	.	+	1
            chr1	6	8	.	.	-	2

        Does not include noncanonical bed-column `Gene`:

        >>> gr.to_bed(keep=False)
        'chr1\t1\t5\t.\t.\t+\nchr1\t6\t8\t.\t.\t-\n'

        File contents::

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
        map_cols: dict | None = None,
    ) -> str | None:
        r"""Write to General Feature Format 3.

        The GFF format consists of a tab-separated file without header.
        GFF contains a fixed amount of columns, indicated below (names before ":").
        For each of these, PyRanges will use the corresponding column (names after ":").

        * seqname: Chromosome
        * source: Source
        * feature: Feature
        * start: Start
        * end: End
        * score: Score
        * strand: Strand
        * phase: Frame
        * attribute: auto-filled

        Columns which are not mapped to GFF columns are appended as a field
        in the ``attribute`` string (i.e. the last field).

        Parameters
        ----------
        path : str, default None, i.e. return string representation.
            Where to write file.

        compression : {'infer', 'gzip', 'bz2', 'zip', 'xz', None}, default "infer"
            Which compression to use. Uses file extension to infer by default.

        map_cols: dict, default None
            Override mapping between GTF and PyRanges fields for any number of columns.
            Format: ``{gtf_column : pyranges_column}``
            If a mapping is found for the "attribute"` column, it is not auto-filled


        Note
        ----
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

        How the file would look::

            1	.	gene	2	4	.	.	.
            1	.	exon	4	6	.	.	.
            1	.	exon	6	9	.	.	.

        >>> gr["Gene"] = [1, 2, 3]
        >>> gr["function"] = ["a b", "c", "def"]
        >>> gr.to_gff3()
        '1\t.\tgene\t2\t4\t.\t.\t.\tGene=1;function=a b\n1\t.\texon\t4\t6\t.\t.\t.\tGene=2;function=c\n1\t.\texon\t6\t9\t.\t.\t.\tGene=3;function=def\n'

        How the file would look::

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

        How the file would look::

            1	.	mRNA	2	4	.	.	0	Gene=1;function=a b
            1	.	CDS	4	6	.	.	2	Gene=2;function=c
            1	.	CDS	6	9	.	.	1	Gene=3;function=def

        >>> gr['custom'] = ['AA', 'BB', 'CC']
        >>> gr
          index  |      Chromosome    Start      End  Feature       Gene  function      phase  custom
          int64  |           int64    int64    int64  object       int64  object        int64  object
        -------  ---  ------------  -------  -------  ---------  -------  ----------  -------  --------
              0  |               1        1        4  mRNA             1  a b               0  AA
              1  |               1        3        6  CDS              2  c                 2  BB
              2  |               1        5        9  CDS              3  def               1  CC
        PyRanges with 3 rows, 8 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> print(gr.to_gff3(map_cols={"feature": "custom"})) # doctest: +NORMALIZE_WHITESPACE
        1	.	AA	2	4	.	.	0	Feature=mRNA;Gene=1;function=a b
        1	.	BB	4	6	.	.	2	Feature=CDS;Gene=2;function=c
        1	.	CC	6	9	.	.	1	Feature=CDS;Gene=3;function=def
        <BLANKLINE>

        >>> print(gr.to_gff3(map_cols={"attribute": "custom"})) # doctest: +NORMALIZE_WHITESPACE
        1	.	mRNA	2	4	.	.	0	AA
        1	.	CDS	4	6	.	.	2	BB
        1	.	CDS	6	9	.	.	1	CC
        <BLANKLINE>

        """
        from pyranges.core.out import _to_gff_like

        return _to_gff_like(self, path=path, out_format="gff3", compression=compression, map_cols=map_cols)

    def to_gtf(
        self,
        path: None = None,
        compression: PANDAS_COMPRESSION_TYPE = None,
        map_cols: dict | None = None,
    ) -> str | None:
        r"""Write to Gene Transfer Format.

        The GTF format consists of a tab-separated file without header.
        It contains a fixed amount of columns, indicated below (names before ":").
        For each of these, PyRanges will use the corresponding column (names after ":").

        * seqname: Chromosome
        * source: Source
        * feature: Feature
        * start: Start
        * end: End
        * score: Score
        * strand: Strand
        * frame: Frame
        * attribute: auto-filled

        Columns which are not mapped to GTF columns are appended as a field
        in the ``attribute`` string (i.e. the last field).

        Parameters
        ----------
        path : str, default None, i.e. return string representation.
            Where to write file.

        compression : {'infer', 'gzip', 'bz2', 'zip', 'xz', None}, default "infer"
            Which compression to use. Uses file extension to infer by default.

        map_cols: dict, default None
            Override mapping between GTF and PyRanges fields for any number of columns.
            Format: ``{gtf_column : pyranges_column}``
            If a mapping is found for the "attribute"` column, it is not auto-filled

        Note
        ----
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

        What the file contents look like::

            1	.	gene	2	4	.	.	.
            1	.	exon	4	6	.	.	.
            1	.	exon	6	9	.	.	.

        >>> gr.Feature = ["GENE", "EXON", "EXON"]
        >>> gr.to_gtf()  # the raw string output
        '1\t.\tGENE\t2\t4\t.\t.\t.\t\n1\t.\tEXON\t4\t6\t.\t.\t.\t\n1\t.\tEXON\t6\t9\t.\t.\t.\t\n'

        The file would look like::

            1	.	GENE	2	4	.	.	.
            1	.	EXON	4	6	.	.	.
            1	.	EXON	6	9	.	.	.

        >>> gr["tag"] = [11, 22, 33]
        >>> gr
          index  |      Chromosome    Start      End  Feature        tag
          int64  |           int64    int64    int64  object       int64
        -------  ---  ------------  -------  -------  ---------  -------
              0  |               1        1        4  GENE            11
              1  |               1        3        6  EXON            22
              2  |               1        5        9  EXON            33
        PyRanges with 3 rows, 5 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> print(gr.to_gff3()) # doctest: +NORMALIZE_WHITESPACE
        1	.	GENE	2	4	.	.	.	tag=11
        1	.	EXON	4	6	.	.	.	tag=22
        1	.	EXON	6	9	.	.	.	tag=33
        <BLANKLINE>

        >>> print(gr.to_gff3(map_cols={'seqname':'tag'})) # doctest: +NORMALIZE_WHITESPACE
        11	.	GENE	2	4	.	.	.	Chromosome=1
        22	.	EXON	4	6	.	.	.	Chromosome=1
        33	.	EXON	6	9	.	.	.	Chromosome=1
        <BLANKLINE>

        >>> print(gr.to_gff3(map_cols={'attribute':'tag'})) # doctest: +NORMALIZE_WHITESPACE
        1	.	GENE	2	4	.	.	.	11
        1	.	EXON	4	6	.	.	.	22
        1	.	EXON	6	9	.	.	.	33
        <BLANKLINE>

        """
        from pyranges.core.out import _to_gff_like

        return _to_gff_like(self, path=path, out_format="gtf", compression=compression, map_cols=map_cols)

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

        """
        broken_examples = (  # noqa:F841
            """

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
        )

        if strand is None:
            strand = self.strand_valid

        from pyranges.methods.to_rle import _to_rle

        strand = validate_and_convert_use_strand(self, strand)
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

        Negative strand intervals are sliced in descending order by default:

        >>> gs = pr.PyRanges({"Chromosome": [1, 1], "Start": [200, 600], "End": [332, 787], "Strand":['+', '-']})
        >>> gs
          index  |      Chromosome    Start      End  Strand
          int64  |           int64    int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |               1      200      332  +
              1  |               1      600      787  -
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> w = gs.window(100)
        >>> w['lengths'] = w.lengths() # add lengths column to see the length of the windows
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
        import ruranges

        use_strand = validate_and_convert_use_strand(self, use_strand)

        negative_strand = (self[STRAND_COL] == "-").to_numpy() if use_strand else np.zeros(len(self), dtype=bool)
        # assert 0, negative_strands
        idx, starts, ends = ruranges.window(  # type: ignore[attr-defined]
            starts=self[START_COL].to_numpy(),
            ends=self[END_COL].to_numpy(),
            negative_strand=negative_strand,
            window_size=window_size,
        )
        df = self.take(idx)  # type: ignore[arg-type]
        df.loc[:, START_COL] = starts
        df.loc[:, END_COL] = ends

        # every interval can be processed individually. This may be optimized in the future.
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
            Column(s) to return.

        preserve_loc_order : bool, default False
            Whether to preserve the order of the genome location columns.
            If False, the genome location columns will be moved to the left.

        Returns
        -------
        PyRanges

            PyRanges with the requested columns.

        See Also
        --------
        PyRanges.remove_nonloc_columns : remove all columns that are not genome location columns.

        Examples
        --------
        >>> gr = pr.PyRanges({"Chromosome": [1], "Start": [895], "Strand": ["+"],
        ...                   "Score": [1], "Score2": [2], "End": [1259]})
        >>> gr
          index  |      Chromosome    Start  Strand      Score    Score2      End
          int64  |           int64    int64  object      int64     int64    int64
        -------  ---  ------------  -------  --------  -------  --------  -------
              0  |               1      895  +               1         2     1259
        PyRanges with 1 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        Genomic location columns are moved to the left by default:

        >>> gr.get_with_loc_columns(["Score2", "Score", "Score2"])
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
        keys = arg_to_list(key)

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

    def group_cumsum(
        self,
        *,
        match_by: VALID_BY_TYPES = None,
        use_strand: VALID_USE_STRAND_TYPE = USE_STRAND_DEFAULT,
        cumsum_start_column: str | None = None,
        cumsum_end_column: str | None = None,
        sort: bool = True,
    ) -> "PyRanges":
        """Strand-aware cumulative length of every interval *within each chromosome-level group*.

        For every chromosome (and, if supplied, every unique combination in
        *match_by*) the intervals are walked 5→3 **on their own strand**.
        Two new columns are added:

        * ``cumsum_start_column`` - running total **before** the interval
        * ``cumsum_end_column``   - running total **after**  the interval

        Parameters
        ----------
        match_by : str or list, default *None*
            Additional column(s) that must match for two intervals to share a
            cumulative coordinate space.  When *None* all intervals on the same
            chromosome are cumulated together.
        cumsum_start_column, cumsum_end_column : str | None, default None
            Names of the columns added to the returned frame. If None is given,
            Start and End is used.
        use_strand: {"auto", True, False}, default: "auto"
            Whether negative strand intervals should be sliced in descending order, meaning 5' to 3'.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).
        sort : bool, default True
            Whether to sort the results by the original row order.

        Returns
        -------
        PyRanges
            Copy of *self* with the two cumulative-length columns appended.

        Examples
        --------
        >>> gr = pr.example_data.ensembl_gtf.get_with_loc_columns(["Feature", "gene_name"])
        >>> gr = gr[gr.Feature == "exon"]
        >>> gr
          index  |      Chromosome    Start      End  Strand      Feature     gene_name
          int64  |        category    int64    int64  category    category    object
        -------  ---  ------------  -------  -------  ----------  ----------  -----------
              2  |               1    11868    12227  +           exon        DDX11L1
              3  |               1    12612    12721  +           exon        DDX11L1
              4  |               1    13220    14409  +           exon        DDX11L1
              5  |               1   112699   112804  -           exon        AL627309.1
              6  |               1   110952   111357  -           exon        AL627309.1
              8  |               1   133373   133723  -           exon        AL627309.1
              9  |               1   129054   129223  -           exon        AL627309.1
             10  |               1   120873   120932  -           exon        AL627309.1
        PyRanges with 8 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.
        >>> gr.group_cumsum(match_by="gene_name")
          index  |      Chromosome    Start      End  Strand      Feature     gene_name
          int64  |        category    int64    int64  category    category    object
        -------  ---  ------------  -------  -------  ----------  ----------  -----------
              2  |               1        0      359  +           exon        DDX11L1
              3  |               1      359      468  +           exon        DDX11L1
              4  |               1      468     1657  +           exon        DDX11L1
              5  |               1      578      683  -           exon        AL627309.1
              6  |               1      683     1088  -           exon        AL627309.1
              8  |               1        0      350  -           exon        AL627309.1
              9  |               1      350      519  -           exon        AL627309.1
             10  |               1      519      578  -           exon        AL627309.1
        PyRanges with 8 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        """
        import ruranges  # local import reduces start-up time

        strand = validate_and_convert_use_strand(self, use_strand)
        match_by = arg_to_list(match_by)
        group_ids = factorize(self, match_by)

        forward = (self[STRAND_COL] == FORWARD_STRAND).to_numpy() if strand else np.ones(self.shape[0], dtype=np.bool_)

        idx, cumsum_start, cumsum_end = ruranges.group_cumsum(  # type: ignore[attr-defined]
            starts=self[START_COL].to_numpy(),
            ends=self[END_COL].to_numpy(),
            groups=group_ids,
            negative_strand=forward,
            sort=sort,
        )

        res = self.take(idx).copy()  # type: ignore[arg-type]
        if cumsum_start_column is None:
            res.loc[:, START_COL] = cumsum_start
            res.loc[:, END_COL] = cumsum_end
        else:
            res.insert(res.shape[1], cumsum_start_column, cumsum_start)
            res.insert(res.shape[1], cumsum_end_column, cumsum_end)

        return mypy_ensure_pyranges(res)

    def intersect(  # type: ignore[override]
        self,
        other: "PyRanges",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto",
        *,
        multiple: VALID_OVERLAP_TYPE = "all",
        match_by: VALID_BY_TYPES = None,
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
            A PyRanges with overlapping intervals. Input index is preserved, but may contain duplicates.

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
              1  |    chr1               28       30  b
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        from pyranges.methods.overlap import _intersect

        # note: argument multiple = 'containment' is formally accepted but omitted in docstring since the result
        # will be always the same as self.overlap(other, contained=True), no intersect is done in that case

        _other, by = prepare_by_binary(
            self,
            other,
            strand_behavior=strand_behavior,
            match_by=match_by,
        )

        gr = _intersect(
            self,
            _other,
            by=by,
            multiple=multiple,
        )

        return mypy_ensure_pyranges(gr)

    def compute_interval_metrics(
        self,
        metrics: str | Iterable[str] | Mapping[str, str] = "fraction",
        *,
        start: str = START_COL,
        end: str = END_COL,
        start2: str = START_COL + JOIN_SUFFIX,
        end2: str = END_COL + JOIN_SUFFIX,
        denom: str = "first",
    ) -> "pr.PyRanges":
        """Attach interval-relationship metrics as new columns.

        Parameters
        ----------
        metrics
            One of the following forms:
              * single string, eg "length"
              * iterable of strings, eg ["fraction", "jaccard"]
              * mapping {metric_name -> new_column_name} to rename on the fly
            Accepted metric names are listed in VALID_METRICS.
        denom
            Denominator for the *fraction* metric.  Must be "first", "second" or "union".
        start, end : str, default START_COL / END_COL
            Column names holding the first interval coordinates.
        start2, end2 : str, default START_COL + "_b" / END_COL + "_b"
            Column names holding the second interval coordinates.
        denom : {"first", "second", "union"}, default "first"
            Denominator used by the *fraction* metric.

        Returns
        -------
        RangeFrame
            Copy of self with extra metric columns.

        Metrics
        -------
        overlap_length
            Raw number of overlapping bases.

        fraction
            Overlap divided by a denominator chosen with *denom*
            ("first", "second", or "union").

        jaccard
            Overlap divided by the union length of the two intervals.

        distance
            Positive gap in bases when intervals do not touch;
            0 when they overlap or abut.

        overlap
            Boolean flag - True if at least one base overlaps.

        signed_distance
            Same as *distance* but signed:
            negative when the second interval is upstream of the first,
            positive when downstream, 0 when touching/overlapping.

        midpoint_distance
            Absolute distance between interval midpoints.

        symmetric_coverage
            2 * overlap ÷ (length1 + length2).  Ranges from 0 to 1.

        relative_direction
            For frames that contain "Strand" and "Strand_b":
            "same" if strands match, "opposite" if they differ,
            "unknown" if either strand is "." or missing.

        Examples
        --------
        >>> import pyranges as pr
        >>> df = pd.DataFrame(
        ...     {
        ...         "Chromosome": ["chr1"] * 5,
        ...         "Start":      [2, 10, 20, 40, 80],
        ...         "End":        [8, 12, 25, 45, 85],
        ...         "Strand":     ["+", "-", "+", "+", "-"],
        ...         "Start_b":    [5,  9, 23, 60, 70],
        ...         "End_b":      [7, 20, 30, 70, 75],
        ...         "Strand_b":   ["+", "+", "-", "-", "+"],
        ...     }
        ... )
        >>> gr = pr.PyRanges(df)

        # length
        >>> gr.compute_interval_metrics("overlap_length")["overlap_length"].tolist()
        [2, 2, 2, 0, 0]

        # fraction (overlap / first interval length)
        >>> gr.compute_interval_metrics("fraction")["fraction"].round(2).tolist()
        [0.33, 1.0, 0.4, 0.0, 0.0]

        # jaccard
        >>> gr.compute_interval_metrics("jaccard")["jaccard"].round(2).tolist()
        [0.33, 0.18, 0.2, 0.0, 0.0]

        # distance (unsigned gap; 0 when overlapping)
        >>> gr.compute_interval_metrics("distance")["distance"].tolist()
        [0, 0, 0, 15, 5]

        # overlap flag
        >>> gr.compute_interval_metrics("overlap")["overlap"].tolist()
        [True, True, True, False, False]

        # signed_distance
        >>> gr.compute_interval_metrics("signed_distance")["signed_distance"].tolist()
        [0, 0, 0, 15, -5]

        # midpoint_distance
        >>> gr.compute_interval_metrics("midpoint_distance")["midpoint_distance"].tolist()
        [1.0, 3.5, 4.0, 22.5, 10.0]

        # symmetric_coverage
        >>> gr.compute_interval_metrics("symmetric_coverage")["symmetric_coverage"].round(2).tolist()
        [0.5, 0.31, 0.33, 0.0, 0.0]

        # relative_direction (requires strand columns)
        >>> gr.compute_interval_metrics("relative_direction")["relative_direction"].tolist()
        ['same', 'opposite', 'opposite', 'opposite', 'opposite']

        """
        return mypy_ensure_pyranges(
            compute_interval_metrics(
                self,
                metrics=metrics,
                start=start,
                end=end,
                start2=start2,
                end2=end2,
                denom=denom,
            )
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

        The function is designed as post-processing after join_ranges to aggregate the coordinates of the two intervals.
        By default, the new start and end columns will be the intersection of the intervals.

        Parameters
        ----------
        function : {"intersect", "union", "swap"} or Callable, default "intersect"
            How to combine the self and other intervals: "intersect", "union", or "swap"
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
              0  |    chr1             9916    10115  -           chr1                 9988    10187  -
              0  |    chr1             9916    10115  -           chr1                10079    10278  -
              1  |    chr1             9939    10138  +           chr1                10073    10272  +
              2  |    chr1             9951    10150  -           chr1                 9988    10187  -
              2  |    chr1             9951    10150  -           chr1                10079    10278  -
        PyRanges with 5 rows, 8 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes and 2 strands.

        The default operation is to intersect the intervals:

        >>> j.combine_interval_columns()
          index  |    Chromosome      Start      End  Strand      Chromosome_b    Strand_b
          int64  |    category        int64    int64  category    category        category
        -------  ---  ------------  -------  -------  ----------  --------------  ----------
              0  |    chr1             9988    10115  -           chr1            -
              0  |    chr1            10079    10115  -           chr1            -
              1  |    chr1            10073    10138  +           chr1            +
              2  |    chr1             9988    10150  -           chr1            -
              2  |    chr1            10079    10150  -           chr1            -
        PyRanges with 5 rows, 6 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes and 2 strands.

        Take the union instead:

        >>> j.combine_interval_columns('union')
          index  |    Chromosome      Start      End  Strand      Chromosome_b    Strand_b
          int64  |    category        int64    int64  category    category        category
        -------  ---  ------------  -------  -------  ----------  --------------  ----------
              0  |    chr1             9916    10187  -           chr1            -
              0  |    chr1             9916    10278  -           chr1            -
              1  |    chr1             9939    10272  +           chr1            +
              2  |    chr1             9951    10187  -           chr1            -
              2  |    chr1             9951    10278  -           chr1            -
        PyRanges with 5 rows, 6 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes and 2 strands.

        >>> j.combine_interval_columns('swap')
          index  |    Chromosome      Start      End  Strand      Chromosome_b    Strand_b
          int64  |    category        int64    int64  category    category        category
        -------  ---  ------------  -------  -------  ----------  --------------  ----------
              0  |    chr1             9988    10187  -           chr1            -
              0  |    chr1            10079    10278  -           chr1            -
              1  |    chr1            10073    10272  +           chr1            +
              2  |    chr1             9988    10187  -           chr1            -
              2  |    chr1            10079    10278  -           chr1            -
        PyRanges with 5 rows, 6 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes and 2 strands.

        Use a custom function that keeps the start of the first interval and the end of the second:

        >>> def custom_combine(s1, e1, s2, e2): return (s1, e2)
        >>> j.combine_interval_columns(custom_combine)
          index  |    Chromosome      Start      End  Strand      Chromosome_b    Strand_b
          int64  |    category        int64    int64  category    category        category
        -------  ---  ------------  -------  -------  ----------  --------------  ----------
              0  |    chr1             9916    10187  -           chr1            -
              0  |    chr1             9916    10278  -           chr1            -
              1  |    chr1             9939    10272  +           chr1            +
              2  |    chr1             9951    10187  -           chr1            -
              2  |    chr1             9951    10278  -           chr1            -
        PyRanges with 5 rows, 6 columns, and 1 index columns (with 2 index duplicates).
        Contains 1 chromosomes and 2 strands.

        """
        res = super().combine_interval_columns(
            function=function,
            start=start,
            end=end,
            start2=start2,
            end2=end2,
            drop_old_columns=drop_old_columns,
        )
        return mypy_ensure_pyranges(res)

    def complement(
        self: "PyRanges",
        match_by: VALID_BY_TYPES = None,
        *,
        use_strand: VALID_USE_STRAND_TYPE = USE_STRAND_DEFAULT,
        include_first_interval: bool = False,
        group_sizes_col: str = CHROM_COL,
        chromsizes: "dict[str | int, int] | None" = None,
    ) -> "PyRanges":
        """Return the internal complement of the intervals, i.e. its introns.

        The complement of an interval is the set of intervals that are not covered by the original interval.
        This function is useful for obtaining the introns of a set of exons, corresponding to the
        "internal" complement, i.e. excluding the first and last portion of each chromosome not covered by intervals.

        Parameters
        ----------
        match_by : str or list, optional
            Column(s) to group intervals (e.g. exons). If provided, the complement will be calculated separately for each group.
        use_strand : {"auto", True, False}, default "auto"
            Whether to return complement intervals separately for those on the positive and negative strands.
            The default "auto" means that strand information is used if present and valid (see .strand_valid).
        include_first_interval : bool, default False
            If True, include the external complement interval at the beginning of the chromosome (or group),
            i.e. the interval from the start of the chromosome up to the first interval.
        group_sizes_col : str, default CHROM_COL
            The column name used to match keys in the ``chromsizes`` mapping. This determines the total size
            of each chromosome (or group) when calculating external complement intervals.
        chromsizes : dict[str | int, int] or None, optional
            If provided, external complement intervals will also be returned, i.e. the intervals corresponding to the
            beginning of the chromosome up to the first interval and from the last interval to the end of the chromosome.
            The dictionary should map chromosome (or group) identifiers to their total sizes. A PyRanges or pyfaidx.Fasta
            object is also accepted since it conveniently loads chromosome lengths.

        Notes
        -----
        * To ensure non-overlap among the input intervals, merge_overlaps is run before the complement is calculated.
        * Bookended intervals will result in no complement intervals returned since they would be of length 0.

        See Also
        --------
        PyRanges.subtract_ranges : report non-overlapping subintervals
        PyRanges.boundaries : report the boundaries of groups of intervals (e.g. transcripts/genes)


        Examples
        --------
        >>> a = pr.PyRanges(dict(Chromosome="chr1", Start=[2, 10, 20, 40], End=[5, 18, 30, 46], ID=['a', 'a', 'b', 'b']))
        >>> a
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                2        5  a
              1  |    chr1               10       18  a
              2  |    chr1               20       30  b
              3  |    chr1               40       46  b
        PyRanges with 4 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> a.complement('ID', group_sizes_col="ID", chromsizes={"a": 22, "b": 100}, include_first_interval=True)
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                0        2  a
              1  |    chr1                5       10  a
              2  |    chr1               18       22  a
              3  |    chr1                0       20  b
              4  |    chr1               30       40  b
              5  |    chr1               46      100  b
        PyRanges with 6 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        Get complement of the whole set of intervals, without grouping:

        Using complement to get introns:

        >>> a.complement('ID')
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                5       10  a
              1  |    chr1               30       40  b
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> a.complement()
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                5       10
              1  |    chr1               18       20
              2  |    chr1               30       40
        PyRanges with 3 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        Include external intervals:

        >>> a.complement(chromsizes={'chr1': 10000}, include_first_interval=True)
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                0        2
              1  |    chr1                5       10
              2  |    chr1               18       20
              3  |    chr1               30       40
              4  |    chr1               46    10000
        PyRanges with 5 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> a.complement('ID', chromsizes={'chr1': 10000}, include_first_interval=True)
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                0        2  a
              1  |    chr1                5       10  a
              2  |    chr1               18    10000  a
              3  |    chr1                0       20  b
              4  |    chr1               30       40  b
              5  |    chr1               46    10000  b
        PyRanges with 6 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        For complement of whole sets of intervals, you can explicitly use_strand or not:

        >>> b = pr.PyRanges(dict(Chromosome="chr1", Start=[1, 10, 20, 40], End=[5, 18, 30, 46],
        ...                      Strand=['+', '+', '-', '-']))
        >>> b
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        5  +
              1  |    chr1               10       18  +
              2  |    chr1               20       30  -
              3  |    chr1               40       46  -
        PyRanges with 4 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> b.complement(use_strand=True)  # same as b.complement() because b.strand_valid == True
          index  |    Chromosome      Start      End  Strand
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                5       10  +
              1  |    chr1               30       40  -
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> b.complement(use_strand=False)
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                5       10
              1  |    chr1               18       20
              2  |    chr1               30       40
        PyRanges with 3 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> b.complement(use_strand=False, chromsizes={'chr1': 10000}, include_first_interval=True)
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                0        1
              1  |    chr1                5       10
              2  |    chr1               18       20
              3  |    chr1               30       40
              4  |    chr1               46    10000
        PyRanges with 5 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        Bookended intervals (indices 0-1 below) and overlapping intervals (2-3) won't return any in-between intervals:

        >>> c = pr.PyRanges(dict(Chromosome="chr1", Start=[1, 5, 8, 10], End=[5, 7, 14, 16]))
        >>> c.complement()
          index  |    Chromosome      Start      End
          int64  |    object          int64    int64
        -------  ---  ------------  -------  -------
              0  |    chr1                7        8
        PyRanges with 1 rows, 3 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        by = prepare_by_single(self, use_strand=use_strand, match_by=match_by)

        result = _complement(
            self,
            by=by,
            chromsizes=chromsizes,
            chromsizes_col=group_sizes_col,
            include_first_interval=include_first_interval,
        )

        return mypy_ensure_pyranges(result)

    def complement_overlaps(  # type: ignore[override]
        self: "PyRanges",
        other: "PyRanges",
        *,
        match_by: VALID_BY_TYPES = None,
        slack: int = 0,
        use_strand: VALID_USE_STRAND_TYPE = USE_STRAND_DEFAULT,
    ) -> "PyRanges":
        """Return the internal complement of self with respect to another PyRanges.

        The complement of an interval is defined as the subinterval(s) in self that are not overlapped by any
        intervals in other. This function is useful for obtaining regions such as introns by excluding intervals
        (e.g. exons) that overlap with those in another PyRanges object.

        Parameters
        ----------
        other : PyRanges
            A PyRanges object whose intervals are used to determine the overlaps.
        match_by : str or list, default None
            Column(s) to group intervals by (e.g. exons). If provided, the complement will be calculated separately
            for each group.
        slack : int, default 0
            An integer offset that adjusts the overlap threshold when computing the complement intervals.
            Negative values reduce the required gap between intervals (effectively "shrinking" them), while positive values
            increase the gap threshold.
        use_strand : {"auto", True, False}, default "auto"
            Whether to compute the complement separately for intervals on the positive and negative strand.
            The default "auto" means that strand information is used if present and valid (see .strand_valid).

        Returns
        -------
        PyRanges
            A PyRanges object containing the complement intervals (i.e. regions in self that do not overlap with any intervals in other).

        Notes
        -----
        * To ensure non-overlap among the input intervals, merge_overlaps is run before the complement is calculated.
        * Bookended intervals will result in no complement intervals returned since they would be of length 0.

        See Also
        --------
        PyRanges.subtract_ranges : Report non-overlapping subintervals.
        PyRanges.boundaries : Report the boundaries of groups of intervals (e.g. transcripts/genes).

        Examples
        --------
        >>> gr = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [1, 4, 10],
        ...                    "End": [3, 9, 11], "ID": ["a", "b", "c"]})
        >>> gr2 = pr.PyRanges({"Chromosome": ["chr1"] * 3, "Start": [2, 2, 9],
        ...                    "End": [3, 9, 10]})
        >>> gr.complement_overlaps(gr2)
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              2  |    chr1               10       11  c
        PyRanges with 1 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        >>> gr.complement_overlaps(gr2, slack=-1)
          index  |    Chromosome      Start      End  ID
          int64  |    object          int64    int64  object
        -------  ---  ------------  -------  -------  --------
              0  |    chr1                1        3  a
              2  |    chr1               10       11  c
        PyRanges with 2 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes.

        """
        result = super().complement_overlaps(
            other=other,
            match_by=prepare_by_single(self, use_strand=use_strand, match_by=match_by),
            slack=slack,
        )

        return mypy_ensure_pyranges(result)

    def get_sequence(
        self: "PyRanges",
        path: Path | None = None,
        pyfaidx_fasta: Optional["pyfaidx.Fasta"] = None,
        use_strand: VALID_USE_STRAND_TYPE = USE_STRAND_DEFAULT,
        sequence_column: str = "Sequence",
    ) -> pd.Series:
        r"""Get the sequence of the intervals from a fasta file.

        Parameters
        ----------
        path : Path
            Path to fasta file. It will be indexed using pyfaidx if an index is not found

        pyfaidx_fasta : pyfaidx.Fasta
            Alternative method to provide fasta target, as a pyfaidx.Fasta object

        use_strand: {"auto", True, False}, default: "auto"
            If True, intervals on the reverse strand will be reverse complemented.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        sequence_column: str, default "Sequence"
            What the added column will be called.

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
        pyranges.seqs : submodule with sequence-related functions

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

        >>> seq = gr.get_sequence("temp.fasta", sequence_column="Sequence")
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

        >>> gr.get_sequence("temp.fasta", use_strand=False)
        0      CAT
        1    GTAAT
        Name: Sequence, dtype: object

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

        use_strand = validate_and_convert_use_strand(self, use_strand=use_strand)

        iterables = (
            zip(self[CHROM_COL], self[START_COL], self[END_COL], [FORWARD_STRAND] * len(self), strict=True)
            if not use_strand
            else zip(self[CHROM_COL], self[START_COL], self[END_COL], self[STRAND_COL], strict=True)
        )

        # below, -seq from pyfaidx is used to get reverse complement of the sequence
        seqs = []
        for chromosome, start, end, strand in iterables:
            _fasta = pyfaidx_fasta[chromosome]
            forward_strand = strand == FORWARD_STRAND
            if (seq := _fasta[start:end]) is not None:
                seqs.append(seq.seq if forward_strand else (-seq).seq)

        return pd.Series(data=seqs, index=self.index, name=sequence_column)

    def get_transcript_sequence(
        self: "PyRanges",
        transcript_id: str | list[str],
        path: Path | None = None,
        *,
        pyfaidx_fasta: Optional["pyfaidx.Fasta"] = None,
        use_strand: VALID_USE_STRAND_TYPE = USE_STRAND_DEFAULT,
        sequence_column: str = "Sequence",
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

        use_strand: {"auto", True, False}, default: "auto"
            If True, intervals on the reverse strand will be reverse complemented.
            The default "auto" means True if PyRanges has valid strands (see .strand_valid).

        sequence_column: str, default "Sequence"
            What the added column will be called.

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
        pyranges.seqs : submodule with sequence-related functions

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

        With use_strand=False, all intervals are treated as if on the forward strand:

        >>> seq2 = gr.get_transcript_sequence(path="temp.fasta", transcript_id='transcript', use_strand=False, sequence_column="Seq2")
        >>> seq2
          transcript     Seq2
        0         t1     AAAC
        1         t2  GGGATTT
        2         t4     GGGA
        3         t5      TTT

        To write to a file in fasta format:
        >>> with open('outfile.fasta', 'w') as fw:
        ...     nchars=60
        ...     for row in seq.itertuples():
        ...         s = '\\n'.join([ row.Sequence[i:i+nchars] for i in range(0, len(row.Sequence), nchars)])
        ...         _bytes_written = fw.write(f'>{row.transcript}\\n{s}\\n')

        """
        use_strand = validate_and_convert_use_strand(self, use_strand=use_strand)
        gr = self.sort_ranges(use_strand=use_strand)

        seq = gr.get_sequence(path=path, pyfaidx_fasta=pyfaidx_fasta, use_strand=use_strand)
        gr[sequence_column] = seq.to_numpy()

        return gr.groupby(transcript_id, as_index=False).agg({sequence_column: "".join})

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
        ValueError: Not all chromosomes were in the chromsize dict.
        Missing keys: {3}.

        """
        import ruranges

        if isinstance(chromsizes, pd.DataFrame):
            chromsizes = dict(zip(chromsizes[CHROM_COL], chromsizes[END_COL], strict=True))
        elif not isinstance(chromsizes, dict):
            # fall-back for pyfaidx.Fasta etc.
            faidx = cast("dict[str | int, list]", chromsizes)
            chromsizes = {k: len(faidx[k]) for k in faidx}

        if missing := set(self[CHROM_COL]) - chromsizes.keys():
            msg = f"Not all chromosomes were in the chromsize dict.\nMissing keys: {missing}."
            raise ValueError(msg)

        chrom_series: pd.Series = self[CHROM_COL]

        codes, uniques = pd.factorize(chrom_series, sort=False)
        codes = codes.astype(np.uint32)

        lengths_per_code = np.fromiter(
            (chromsizes[u] for u in uniques),
            dtype=np.int64,
            count=len(uniques),
        )

        chrom_lengths_vec = lengths_per_code[codes]

        idxs, starts, ends = ruranges.genome_bounds(
            groups=codes,  # type: ignore[arg-type]
            starts=self[START_COL].to_numpy(),
            ends=self[END_COL].to_numpy(),
            chrom_length=chrom_lengths_vec,  # **row-aligned!**
            clip=clip,
            only_right=only_right,
        )

        res = self.take(idxs)  # type: ignore[arg-type]
        if clip:
            res.loc[:, START_COL] = starts
            res.loc[:, END_COL] = ends

        return mypy_ensure_pyranges(res)


def process_site(df: "PyRanges", slack: int = 0, site: str = "tss") -> "PyRanges":
    """Process a DataFrame to update the 'Start' and 'End' positions based on the site type.

    Args:
    ----
        df: A DataFrame that contains 'Start', 'End', and 'Strand' columns.
        slack: An integer value to adjust the positions.
        site: A string indicating the site type: "tss" or "tes".
              For "tss": use 'Start' for '+' strand and 'End - 1' for '-' strand.
              For "tes": use 'End - 1' for '+' strand and 'Start' for '-' strand.

    Returns:
    -------
        A new DataFrame with updated 'Start' and 'End' columns.

    """
    df = df.copy(deep=True)
    dtype = df.dtypes["Start"]

    if site == "tss":
        starts = np.where(df.Strand == "+", df.Start, df.End - 1)
    elif site == "tes":
        starts = np.where(df.Strand == "+", df.End - 1, df.Start)
    else:
        msg = "Invalid site type. Expected 'tss' or 'tes'."
        raise ValueError(msg)

    ends = starts + slack + 1
    starts = starts - slack
    starts = np.where(starts < 0, 0, starts)

    df.loc[:, "Start"] = starts.astype(dtype)
    df.loc[:, "End"] = ends.astype(dtype)
    return df
