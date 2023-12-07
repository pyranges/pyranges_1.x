"""Statistics useful for genomics."""

from collections import defaultdict
from math import sqrt
from typing import Any, Dict, List, Optional, Union, TYPE_CHECKING

import numpy as np
import pandas as pd
from numpy import ndarray
from pandas.core.frame import DataFrame
from pandas.core.series import Series

import pyranges as pr
from pyranges.methods.statistics import _relative_distance
from pyranges.names import VALID_STRAND_BEHAVIOR_TYPE, STRAND_BEHAVIOR_AUTO, STRAND_BEHAVIOR_IGNORE, \
    STRAND_BEHAVIOR_SAME


if TYPE_CHECKING:
    from pyranges import PyRanges


__all__ = [
    "simes",
    "fisher_exact",
    "StatisticsMethods",
    "fdr",
    "rowbased_rankdata",
    "rowbased_pearson",
    "rowbased_spearman",
    "mcc",
]


def fdr(p_vals: Series) -> Series:
    """Adjust p-values with Benjamini-Hochberg.

    Parameters
    ----------
    data : array-like


    Returns
    -------
    Pandas.DataFrame

        DataFrame where values are order of data.

    Examples
    --------
    >>> import pyranges as pr
    >>> d = {'Chromosome': ['chr3', 'chr6', 'chr13'], 'Start': [146419383, 39800100, 24537618], 'End': [146419483, 39800200, 24537718], 'Strand': ['-', '+', '-'], 'PValue': [0.0039591368855297175, 0.0037600512992788937, 0.0075061166500909205]}
    >>> gr = pr.PyRanges(d)
    >>> gr
    Chromosome        Start        End  Strand        PValue
    object            int64      int64  object       float64
    ------------  ---------  ---------  --------  ----------
    chr3          146419383  146419483  -         0.00395914
    chr6           39800100   39800200  +         0.00376005
    chr13          24537618   24537718  -         0.00750612
    PyRanges with 3 rows and 5 columns.
    Contains 3 chromosomes and 2 strands.

    >>> gr.col.FDR = pr.stats.fdr(gr.PValue)
    >>> gr
    Chromosome        Start        End  Strand        PValue         FDR
    object            int64      int64  object       float64     float64
    ------------  ---------  ---------  --------  ----------  ----------
    chr3          146419383  146419483  -         0.00395914  0.00593871
    chr6           39800100   39800200  +         0.00376005  0.0112802
    chr13          24537618   24537718  -         0.00750612  0.00750612
    PyRanges with 3 rows and 6 columns.
    Contains 3 chromosomes and 2 strands.
    """

    from scipy.stats import rankdata  # type: ignore

    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr


def fisher_exact(
    tp: Series, fp: Series, fn: Series, tn: Series, pseudocount: int = 0
) -> DataFrame:
    """Fisher's exact for contingency tables.

    Computes the hypotheses two-sided, less and greater at the same time.

    The odds-ratio is

    Parameters
    ----------
    tp : array-like of int

        Top left square of contingency table (true positives).

    fp : array-like of int

        Top right square of contingency table (false positives).

    fn : array-like of int

        Bottom left square of contingency table (false negatives).

    tn : array-like of int

        Bottom right square of contingency table (true negatives).

    pseudocount : float, default 0

        Values > 0 allow Odds Ratio to always be a finite number.

    Notes
    -----

    The odds-ratio is computed thusly:

    ``((tp + pseudocount) / (fp + pseudocount)) / ((fn + pseudocount) / (tn + pseudocount))``

    Returns
    -------
    pandas.DataFrame

        DataFrame with columns OR and P, PLeft and PRight.

    See Also
    --------

    pr.stats.fdr : correct for multiple testing

    Examples
    --------

    >>> d = {"TP": [12, 0], "FP": [5, 12], "TN": [29, 10], "FN": [2, 2]}
    >>> df = pd.DataFrame(d)
    >>> df
       TP  FP  TN  FN
    0  12   5  29   2
    1   0  12  10   2

    >>> pr.stats.fisher_exact(df.TP, df.FP, df.TN, df.FN)
             OR         P     PLeft    PRight
    0  0.165517  0.080269  0.044555  0.994525
    1  0.000000  0.000067  0.000034  1.000000
    """

    try:
        from fisher import pvalue_npy  # type: ignore
    except ImportError:
        import sys

        print(
            "fisher needs to be installed to use fisher exact. pip install fisher or conda install -c bioconda fisher."
        )
        sys.exit(-1)

    _tp = np.array(tp, dtype=np.uint)
    _fp = np.array(fp, dtype=np.uint)
    _fn = np.array(fn, dtype=np.uint)
    _tn = np.array(tn, dtype=np.uint)

    left, right, twosided = pvalue_npy(_tp, _fp, _fn, _tn)

    OR = ((_tp + pseudocount) / (_fp + pseudocount)) / (
        (_fn + pseudocount) / (_tn + pseudocount)
    )

    df = pd.DataFrame({"OR": OR, "P": twosided, "PLeft": left, "PRight": right})

    return df


def mcc(
    grs: List["PyRanges"],
    genome: Optional[Union["PyRanges", pd.DataFrame, Dict[str, int]]] = None,
    labels: Optional[str] = None,
    strand: bool = False,
    verbose: bool = False,
) -> DataFrame:
    """Compute Matthew's correlation coefficient for PyRanges overlaps.

    Parameters
    ----------
    grs : list of PyRanges

        PyRanges to compare.

    genome : DataFrame or dict, default None

        Should contain chromosome sizes. By default, end position of the
        rightmost intervals are used as proxies for the chromosome size, but
        it is recommended to use a genome.

    labels : list of str, default None

        Names to give the PyRanges in the output.

    strand : bool, default False

        Whether to compute correlations per strand.

    verbose : bool, default False

        Warn if some chromosomes are in the genome, but not in the PyRanges.

    Examples
    --------
    >>> grs = [pr.data.aorta, pr.data.aorta, pr.data.aorta2]
    >>> mcc = pr.stats.mcc(grs, labels="abc", genome={"chr1": 2100000})
    >>> mcc
       T  F   TP   FP       TN   FN      MCC
    0  a  a  728    0  2099272    0  1.00000
    1  a  b  728    0  2099272    0  1.00000
    3  a  c  457  485  2098787  271  0.55168
    2  b  a  728    0  2099272    0  1.00000
    5  b  b  728    0  2099272    0  1.00000
    6  b  c  457  485  2098787  271  0.55168
    4  c  a  457  271  2098787  485  0.55168
    7  c  b  457  271  2098787  485  0.55168
    8  c  c  942    0  2099058    0  1.00000

    To create a symmetric matrix (useful for heatmaps of correlations):

    >>> mcc.set_index(["T", "F"]).MCC.unstack().rename_axis(None, axis=0)
    F        a        b        c
    a  1.00000  1.00000  0.55168
    b  1.00000  1.00000  0.55168
    c  0.55168  0.55168  1.00000
    """

    import sys
    from itertools import chain, combinations_with_replacement

    if genome is None:
        genome = defaultdict(int)
        for gr in grs:
            for k, v in gr:
                genome[k] = max(genome[k], v.End.max())

    if not isinstance(genome, dict):
        _genome = genome
        genome_length = int(_genome.End.sum())
    else:
        _genome = pd.DataFrame(
            {
                "Chromosome": list(genome.keys()),
                "Start": 0,
                "End": list(genome.values()),
            }
        )
        genome_length = sum(genome.values())

    if labels is None:
        _labels = combinations_with_replacement(np.arange(len(grs)), r=2)
    else:
        assert len(labels) == len(grs)
        _labels = combinations_with_replacement(labels, r=2)

    # remove all non-loc columns before computation
    grs = [gr.merge_overlaps(strand=strand) for gr in grs]

    if _genome is not None:
        genome_length = int(_genome.End.sum())

        if verbose:
            # check that genome definition does not have many more
            # chromosomes than datafiles
            gr_cs = set(chain(*[gr.Chromosome for gr in grs]))

            g_cs = set(_genome.keys())
            surplus = g_cs - gr_cs
            if len(surplus):
                print(
                    "The following chromosomes are in the genome, but not the PyRanges:",
                    ", ".join(surplus),
                    file=sys.stderr,
                )

        if strand:

            def make_stranded(df):
                df = df.copy()
                df2 = df.copy()
                df.insert(df.shape[1], "Strand", "+")
                df2.insert(df2.shape[1], "Strand", "-")
                return pd.concat([df, df2])

            _genome = _genome.apply(make_stranded)

    strand_behavior = STRAND_BEHAVIOR_SAME if strand else STRAND_BEHAVIOR_IGNORE

    rowdicts = []
    for (lt, lf), (t, f) in zip(_labels, combinations_with_replacement(grs, r=2)):
        if verbose:
            print(lt, lf, file=sys.stderr)

        if lt == lf:
            if not strand:
                tp = t.length
                fn = 0
                tn = genome_length - tp
                fp = 0
                rowdicts.append(
                    {"T": lt, "F": lf, "TP": tp, "FP": fp, "TN": tn, "FN": fn, "MCC": 1}
                )
            else:
                for _strand in "+ -".split():
                    tp = t[strand].length
                    fn = 0
                    tn = genome_length - tp
                    fp = 0
                    rowdicts.append(
                        {
                            "T": lt,
                            "F": lf,
                            "Strand": _strand,
                            "TP": tp,
                            "FP": fp,
                            "TN": tn,
                            "FN": fn,
                            "MCC": 1,
                        }
                    )
            continue

        else:
            j = t.interval_join(f, strand_behavior=strand_behavior)
            tp_gr = j.intersect_interval_columns(start2="Start_b", end2="End_b").merge_overlaps(strand=strand)
            if strand:
                for _strand in "+ -".split():
                    tp = tp_gr[_strand].length
                    fp = f[_strand].length - tp
                    fn = t[_strand].length - tp
                    tn = genome_length - (tp + fp + fn)
                    mcc = _mcc(tp, fp, tn, fn)
                    rowdicts.append(
                        {
                            "T": lt,
                            "F": lf,
                            "Strand": _strand,
                            "TP": tp,
                            "FP": fp,
                            "TN": tn,
                            "FN": fn,
                            "MCC": mcc,
                        }
                    )
                    rowdicts.append(
                        {
                            "T": lf,
                            "F": lt,
                            "Strand": _strand,
                            "TP": tp,
                            "FP": fn,
                            "TN": tn,
                            "FN": fp,
                            "MCC": mcc,
                        }
                    )
            else:
                tp = tp_gr.length
                fp = f.length - tp
                fn = t.length - tp
                tn = genome_length - (tp + fp + fn)
                mcc = _mcc(tp, fp, tn, fn)

                rowdicts.append(
                    {
                        "T": lt,
                        "F": lf,
                        "TP": tp,
                        "FP": fp,
                        "TN": tn,
                        "FN": fn,
                        "MCC": mcc,
                    }
                )
                rowdicts.append(
                    {
                        "T": lf,
                        "F": lt,
                        "TP": tp,
                        "FP": fn,
                        "TN": tn,
                        "FN": fp,
                        "MCC": mcc,
                    }
                )

    df = pd.DataFrame.from_records(rowdicts).sort_values(["T", "F"])

    return df


def rowbased_spearman(x: ndarray, y: ndarray) -> ndarray:
    """Fast row-based Spearman's correlation.

    Parameters
    ----------
    x : matrix-like

        2D numerical matrix. Same size as y.

    y : matrix-like

        2D numerical matrix. Same size as x.

    Returns
    -------
    numpy.array

        Array with same length as input, where values are P-values.

    See Also
    --------

    pyranges.statistics.rowbased_pearson : fast row-based Pearson's correlation.
    pr.stats.fdr : correct for multiple testing

    Examples
    --------

    >>> x = np.array([[7, 2, 9], [3, 6, 0], [0, 6, 3]])
    >>> y = np.array([[5, 3, 2], [9, 6, 0], [7, 3, 5]])

    Perform Spearman's correlation pairwise on each row in 10x10 matrixes:

    >>> pr.stats.rowbased_spearman(x, y)
    array([-0.5,  0.5, -1. ])
    """

    x = np.asarray(x)
    y = np.asarray(y)

    rx = rowbased_rankdata(x)
    ry = rowbased_rankdata(y)

    return rowbased_pearson(rx, ry)


def rowbased_pearson(
    x: Union[ndarray, DataFrame], y: Union[ndarray, DataFrame]
) -> ndarray:
    """Fast row-based Pearson's correlation.

    Parameters
    ----------
    x : matrix-like

        2D numerical matrix. Same size as y.

    y : matrix-like

        2D numerical matrix. Same size as x.

    Returns
    -------
    numpy.array

        Array with same length as input, where values are P-values.

    See Also
    --------

    pyranges.statistics.rowbased_spearman : fast row-based Spearman's correlation.
    pr.stats.fdr : correct for multiple testing

    Examples
    --------

    >>> x = np.array([[7, 2, 9], [3, 6, 0], [0, 6, 3]])
    >>> y = np.array([[5, 3, 2], [9, 6, 0], [7, 3, 5]])

    Perform Pearson's correlation pairwise on each row in 10x10 matrixes:

    >>> pr.stats.rowbased_pearson(x, y)
    array([-0.09078413,  0.65465367, -1.        ])
    """

    # Thanks to https://github.com/dengemann

    def ss(a, axis):
        return np.sum(a * a, axis=axis)

    x = np.asarray(x)
    y = np.asarray(y)

    mx = x.mean(axis=-1)
    my = y.mean(axis=-1)

    xm, ym = x - mx[..., None], y - my[..., None]

    r_num = np.add.reduce(xm * ym, axis=-1)
    r_den = np.sqrt(ss(xm, axis=-1) * ss(ym, axis=-1))

    with np.errstate(divide="ignore", invalid="ignore"):
        r = r_num / r_den

    return r


def rowbased_rankdata(data: ndarray) -> DataFrame:
    """Rank order of entries in each row.

    Same as SciPy rankdata with method=mean.

    Parameters
    ----------
    data : matrix-like

        The data to find order of.

    Returns
    -------
    Pandas.DataFrame

        DataFrame where values are order of data.

    Examples
    --------

    >>> x = np.random.randint(10, size=(3, 4))
    >>> x = np.array([[3, 7, 6, 0], [1, 3, 8, 9], [5, 9, 3, 5]])
    >>> pr.stats.rowbased_rankdata(x)
         0    1    2    3
    0  2.0  4.0  3.0  1.0
    1  1.0  2.0  3.0  4.0
    2  2.5  4.0  1.0  2.5
    """

    dc = np.asarray(data).copy()
    sorter = np.apply_along_axis(np.argsort, 1, data)

    inv = np.empty(data.shape, np.intp)

    ranks = np.tile(np.arange(data.shape[1]), (len(data), 1))

    np.put_along_axis(inv, sorter, ranks, axis=1)

    dc = np.take_along_axis(dc, sorter, 1)

    res = np.apply_along_axis(lambda r: r[1:] != r[:-1], 1, dc)

    obs = np.column_stack([np.ones(len(res), dtype=bool), res])

    dense = pd.DataFrame(
        np.take_along_axis(np.apply_along_axis(np.cumsum, 1, obs), inv, 1)
    )

    len_r = obs.shape[1]

    nonzero = np.count_nonzero(obs, axis=1)

    _ranks = []
    for _nonzero, nzdf in pd.DataFrame(obs).groupby(pd.Series(nonzero), sort=False):
        nz = np.apply_along_axis(lambda r: np.nonzero(r)[0], 1, nzdf)

        _count = np.column_stack([nz, np.ones(len(nz)) * len_r])
        _dense = dense.reindex(nzdf.index).values

        _result = 0.5 * (
            np.take_along_axis(_count, _dense, 1)
            + np.take_along_axis(_count, _dense - 1, 1)
            + 1
        )

        result = pd.DataFrame(_result, index=nzdf.index)
        _ranks.append(result)

    final = pd.concat(_ranks).sort_index(kind="mergesort")

    return final


def simes(df, groupby, pcol, keep_position=False):
    """Apply Simes method for giving dependent events a p-value.

    Parameters
    ----------
    df : pandas.DataFrame

        Data to analyse with Simes.

    groupby : str or list of str

        Features equal in these columns will be merged with Simes.

    pcol : str

        Name of column with p-values.

    keep_position : bool, default False

        Keep columns "Chromosome", "Start", "End" and "Strand" if they exist.

    See Also
    --------

    pr.stats.fdr : correct for multiple testing

    Examples
    --------

    >>> s = '''Chromosome Start End Strand Gene PValue
    ... 1 10 20 + P53 0.0001
    ... 1 20 35 + P53 0.0002
    ... 1 30 40 + P53 0.0003
    ... 2 60 65 - FOX 0.05
    ... 2 70 75 - FOX 0.0000001
    ... 2 80 90 - FOX 0.0000021'''

    >>> gr = pr.from_string(s)
    >>> gr
      Chromosome    Start      End  Strand    Gene         PValue
           int64    int64    int64  object    object      float64
    ------------  -------  -------  --------  --------  ---------
               1       10       20  +         P53         0.0001
               1       20       35  +         P53         0.0002
               1       30       40  +         P53         0.0003
               2       60       65  -         FOX         0.05
               2       70       75  -         FOX         1e-07
               2       80       90  -         FOX         2.1e-06
    PyRanges with 6 rows and 6 columns.
    Contains 2 chromosomes and 2 strands.

    >>> simes = pr.stats.simes(gr, "Gene", "PValue")
    >>> simes
      Gene         Simes
    0  FOX  3.000000e-07
    1  P53  3.000000e-04

    >>> pr.stats.simes(gr, "Gene", "PValue", keep_position=True)
      Chromosome    Start      End      Simes  Strand    Gene
           int64    int64    int64    float64  object    object
    ------------  -------  -------  ---------  --------  --------
               2       60       90     1e-07   -         FOX
               1       10       40     0.0001  +         P53
    PyRanges with 2 rows and 6 columns.
    Contains 2 chromosomes and 2 strands.
    """

    if isinstance(groupby, str):
        groupby = [groupby]

    positions = []
    if "Strand" in df:
        stranded = True

    if keep_position:
        positions += ["Chromosome", "Start", "End"]
        if stranded:
            positions += ["Strand"]

    sorter = groupby + [pcol]

    sdf = df[positions + sorter].sort_values(sorter)
    g = sdf.groupby(positions + groupby)

    ranks = g.cumcount().values + 1
    size = g.size().values
    size = np.repeat(size, size)
    multiplied = sdf[pcol].values * size

    simes = multiplied / ranks

    sdf.insert(sdf.shape[1], "Simes", simes)

    if keep_position:
        grpby_dict = {
            "Chromosome": "first",
            "Start": "min",
            "End": "max",
            "Simes": "min",
        }

        if stranded:
            grpby_dict["Strand"] = "first"

        simes = sdf.groupby(groupby).agg(grpby_dict).reset_index()
        columns = list(simes.columns)
        columns.append(columns[0])
        del columns[0]
        simes = pr.PyRanges(simes[columns])
    else:
        simes = sdf.groupby(groupby).Simes.min().reset_index()

    return simes


def chromsizes_as_int(chromsizes: "PyRanges | DataFrame | dict[Any, int]") -> int:
    if isinstance(chromsizes, dict):
        _chromsizes = sum(chromsizes.values())
    elif isinstance(chromsizes, (pd.DataFrame, pr.PyRanges)):
        _chromsizes = chromsizes.End.sum()
    else:
        raise TypeError(
            "chromsizes must be dict, DataFrame or PyRanges, was {}".format(
                type(chromsizes)
            )
        )

    return _chromsizes


class StatisticsMethods:

    """Namespace for statistical comparsion-operations.

    Accessed with gr.stats."""

    def __init__(self, pr: "PyRanges") -> None:
        self.pr = pr

    def forbes(
        self,
        other: "PyRanges",
        chromsizes: "PyRanges | DataFrame | dict[Any, int]",
        strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto"
    ) -> float:
        """Compute Forbes coefficient.

        Ratio which represents observed versus expected co-occurence.

        Described in ``Forbes SA (1907): On the local distribution of certain Illinois fishes: an essay in statistical ecology.``

        Parameters
        ----------
        other : PyRanges

            Intervals to compare with.

        chromsizes : int, dict, DataFrame or PyRanges

            Integer representing genome length or mapping from chromosomes
            to its length.

        strandedness : {None, "same", "opposite", False}, default None, i.e. "auto"

            Whether to compute without regards to strand or on same or opposite.

        Returns
        -------
        float

            Ratio of observed versus expected co-occurence.

        See Also
        --------

        pyranges.statistics.jaccard : compute the jaccard coefficient

        Examples
        --------
        >>> gr, gr2 = pr.data.f1, pr.data.f2
        >>> gr.stats.forbes(gr2, chromsizes={"chr1": 10})
        1.6666666666666667
        """

        _chromsizes = chromsizes_as_int(chromsizes)

        self.pr.ensure_strand_behavior_options_valid(other, strand_behavior=strand_behavior)
        strand = self.pr.strand_values_valid and other.strand_values_valid and strand_behavior in [STRAND_BEHAVIOR_AUTO, True]
        reference_length = self.pr.merge_overlaps(strand=strand).length
        query_length = other.merge_overlaps(strand=strand).length

        intersection_sum = (
            self.pr.set_intersect(other, strand_behavior=strand_behavior).lengths().sum()
        )
        forbes = _chromsizes * intersection_sum / (reference_length * query_length)

        return forbes

    def jaccard(self, other: "PyRanges", chromsizes: dict[str, int], strand_behavior: VALID_STRAND_BEHAVIOR_TYPE = "auto") -> float:
        """Compute Jaccards coefficient.

        Ratio of the intersection and union of two sets.

        Parameters
        ----------
        other : PyRanges

            Intervals to compare with.

        chromsizes : int, dict, DataFrame or PyRanges

            Integer representing genome length or mapping from chromosomes
            to its length.

        strandedness : {None, "same", "opposite", False}, default None, i.e. "auto"

            Whether to compute without regards to strand or on same or opposite.

        Returns
        -------
        float

            Ratio of the intersection and union of two sets.

        See Also
        --------

        pyranges.statistics.forbes : compute the forbes coefficient

        Examples
        --------

        >>> gr, gr2 = pr.data.f1, pr.data.f2
        >>> chromsizes = pr.data.chromsizes
        >>> gr.stats.jaccard(gr2, chromsizes=chromsizes)
        0.3333333333333333
        """

        self.pr.ensure_strand_behavior_options_valid(other, strand_behavior=strand_behavior)
        strand = self.pr.strand_values_valid and other.strand_values_valid and strand_behavior in [STRAND_BEHAVIOR_AUTO, True]

        intersection_sum = self.pr.set_intersect(other).lengths().sum()

        union_sum = 0
        for gr in [self.pr, other]:
            union_sum += gr.merge_overlaps(strand=strand).lengths().sum()

        denominator = union_sum - intersection_sum
        if denominator == 0:
            return 1.0
        else:
            jc = intersection_sum / denominator

        return jc

    def relative_distance(self, other: "PyRanges", **kwargs) -> DataFrame:
        """Compute spatial correllation between two sets.

        Metric which describes relative distance between each interval in one
        set and two closest intervals in another.

        Parameters
        ----------
        other : PyRanges

            Intervals to compare with.

        chromsizes : int, dict, DataFrame or PyRanges

            Integer representing genome length or mapping from chromosomes
            to its length.

        strandedness : {None, "same", "opposite", False}, default None, i.e. "auto"

            Whether to compute without regards to strand or on same or opposite.

        Returns
        -------
        pandas.DataFrame

            DataFrame containing the frequency of each relative distance.

        See Also
        --------

        pyranges.statistics.jaccard : compute the jaccard coefficient
        pyranges.statistics.forbes : compute the forbes coefficient

        Examples
        --------

        >>> gr1, gr2 = pr.data.chipseq, pr.data.chipseq_background
        >>> gr = pd.concat([gr1, gr1.head(4), gr2.tail(4)])
        >>> chromsizes = pr.data.chromsizes
        >>> gr.stats.relative_distance(gr2)
            reldist  count  total  fraction
        0      0.00      4     18  0.222222
        1      0.03      1     18  0.055556
        2      0.04      1     18  0.055556
        3      0.10      1     18  0.055556
        4      0.12      1     18  0.055556
        5      0.13      1     18  0.055556
        6      0.19      1     18  0.055556
        7      0.23      1     18  0.055556
        8      0.24      1     18  0.055556
        9      0.38      1     18  0.055556
        10     0.41      2     18  0.111111
        11     0.42      1     18  0.055556
        12     0.43      2     18  0.111111
        """

        result = pd.Series(_relative_distance(self.pr, other))

        not_nan = ~np.isnan(result)
        result.loc[not_nan] = np.floor(result[not_nan] * 100) / 100
        vc = result.value_counts(dropna=False).to_frame().reset_index()
        vc.columns = "reldist count".split()
        vc.insert(vc.shape[1], "total", len(result))
        vc.insert(vc.shape[1], "fraction", vc["count"] / len(result))
        vc = vc.sort_values("reldist", ascending=True)
        vc = vc.reset_index(drop=True)

        return vc


def _mcc(tp: int, fp: int, tn: int, fn: int) -> float:
    # https://stackoverflow.com/a/56875660/992687
    x = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    return ((tp * tn) - (fp * fn)) / sqrt(x)
