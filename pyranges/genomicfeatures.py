import typing
from typing import TYPE_CHECKING, Any

import pandas as pd
from pandas.core.frame import DataFrame

from pyranges.core.names import CHROM_COL, END_COL, START_COL

__all__ = ["genome_bounds", "tile_genome", "GenomicFeaturesMethods"]

from pyranges.core.pyranges_helpers import mypy_ensure_pyranges

if TYPE_CHECKING:
    try:
        import pyfaidx  # type: ignore[import]

        FastaIdx = pyfaidx.Fasta
    except ImportError:
        FastaIdx = Any
    from pyranges import PyRanges


class GenomicFeaturesMethods:
    """Namespace for methods using feature information.

    Accessed through `gr.features`.
    """

    def __init__(self, pr: "PyRanges") -> None:
        self.pr = pr

    def tss(self) -> "PyRanges":
        """Return the transcription start sites.

        Returns the 5' for every interval with feature "transcript".

        See Also
        --------
        pyranges.genomicfeatures.GenomicFeaturesMethods.tes : return the transcription end sites

        Examples
        --------
        >>> import pyranges as pr
        >>> gr = pr.example_data.ensembl_gtf.get_with_loc_columns(["Source", "Feature"])
        >>> gr
        index    |    Chromosome    Start    End      Strand      Source    Feature
        int64    |    category      int64    int64    category    object    category
        -------  ---  ------------  -------  -------  ----------  --------  ----------
        0        |    1             11868    14409    +           havana    gene
        1        |    1             11868    14409    +           havana    transcript
        2        |    1             11868    12227    +           havana    exon
        3        |    1             12612    12721    +           havana    exon
        ...      |    ...           ...      ...      ...         ...       ...
        7        |    1             120724   133723   -           ensembl   transcript
        8        |    1             133373   133723   -           ensembl   exon
        9        |    1             129054   129223   -           ensembl   exon
        10       |    1             120873   120932   -           ensembl   exon
        PyRanges with 11 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.features.tss()
          index  |      Chromosome    Start      End  Strand      Source    Feature
          int64  |        category    int64    int64  category    object    object
        -------  ---  ------------  -------  -------  ----------  --------  ---------
              0  |               1    11868    11869  +           havana    tss
              1  |               1   133722   133723  -           ensembl   tss
        PyRanges with 2 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        """
        gr = self.pr

        if not gr.strand_valid:
            msg = (
                "Cannot compute TSSes or TESes without strand info. Perhaps use extend()"
                "or subsequence() or spliced_subsequence() instead?"
            )
            raise AssertionError(msg)

        gr = mypy_ensure_pyranges(gr.loc[gr.Feature == "transcript"])
        gr = gr.groupby(CHROM_COL).apply(_tss).reset_index(drop=True)

        gr.Feature = "tss"

        return gr

    def tes(self) -> "PyRanges":
        """Return the transcription end sites.

        Returns the 3' for every interval with feature "transcript".

        See Also
        --------
        pyranges.genomicfeatures.GenomicFeaturesMethods.tss : return the transcription start sites

        Examples
        --------
        >>> import pyranges as pr
        >>> gr = pr.example_data.ensembl_gtf.get_with_loc_columns(["Source", "Feature"])
        >>> gr
        index    |    Chromosome    Start    End      Strand      Source    Feature
        int64    |    category      int64    int64    category    object    category
        -------  ---  ------------  -------  -------  ----------  --------  ----------
        0        |    1             11868    14409    +           havana    gene
        1        |    1             11868    14409    +           havana    transcript
        2        |    1             11868    12227    +           havana    exon
        3        |    1             12612    12721    +           havana    exon
        ...      |    ...           ...      ...      ...         ...       ...
        7        |    1             120724   133723   -           ensembl   transcript
        8        |    1             133373   133723   -           ensembl   exon
        9        |    1             129054   129223   -           ensembl   exon
        10       |    1             120873   120932   -           ensembl   exon
        PyRanges with 11 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.features.tes()
          index  |      Chromosome    Start      End  Strand      Source    Feature
          int64  |        category    int64    int64  category    object    object
        -------  ---  ------------  -------  -------  ----------  --------  ---------
              0  |               1    14408    14409  +           havana    tes
              1  |               1   120724   120725  -           ensembl   tes
        PyRanges with 2 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        """
        gr = self.pr

        if not gr.strand_valid:
            msg = "Cannot compute TSSes or TESes without strand info. Perhaps use extend() or subsequence() or spliced_subsequence() instead?"
            raise ValueError(msg)

        _gr = gr[gr.Feature == "transcript"]
        _gr = _gr.groupby(CHROM_COL).apply(_tes).reset_index(drop=True)

        _gr.Feature = "tes"

        return mypy_ensure_pyranges(_gr)

    def introns(
        self,
        feature_column: str = "Feature",
        outer_feature: str = "gene",
        inner_feature: str = "exon",
        by: str | list[str] | None = None,
    ) -> "PyRanges":
        """Return the introns.

        Parameters
        ----------
        feature_column: str, default "Feature"
            Column to use for feature information.

        outer_feature: str, default "gene"
            Feature to use as outer feature (typically transcript in the case of actual exons).

        inner_feature: str, default "exon"
            Feature to use as inner feature (typically exon in the case of actual exons).

        by : str, {"gene", "transcript"}, default "gene"
            Whether to find introns per gene or transcript.

        See Also
        --------
        pyranges.genomicfeatures.GenomicFeaturesMethods.tss : return the transcription start sites

        Examples
        --------
        >>> import pyranges as pr
        >>> gr = pr.example_data.ensembl_gtf.get_with_loc_columns(["Feature", "gene_id", "transcript_id"])
        >>> gr = gr[gr["gene_id"] == "ENSG00000223972"]
        >>> gr
          index  |      Chromosome    Start      End  Strand      Feature     gene_id          transcript_id
          int64  |        category    int64    int64  category    category    object           object
        -------  ---  ------------  -------  -------  ----------  ----------  ---------------  ---------------
              0  |               1    11868    14409  +           gene        ENSG00000223972  nan
              1  |               1    11868    14409  +           transcript  ENSG00000223972  ENST00000456328
              2  |               1    11868    12227  +           exon        ENSG00000223972  ENST00000456328
              3  |               1    12612    12721  +           exon        ENSG00000223972  ENST00000456328
              4  |               1    13220    14409  +           exon        ENSG00000223972  ENST00000456328
        PyRanges with 5 rows, 7 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        >>> gr = pr.from_string('''Chromosome   Start      End  Strand      Feature Id
        ...          1   0    100  +           gene A
        ...          1   10    20  +           exon A
        ...          1   35    45  +           exon A
        ...          1   30    40  +           exon A
        ...          1   0    50   +           gene B
        ...          1   20   30  +           exon  B''')
        >>> gr
          index  |      Chromosome    Start      End  Strand    Feature    Id
          int64  |           int64    int64    int64  object    object     object
        -------  ---  ------------  -------  -------  --------  ---------  --------
              0  |               1        0      100  +         gene       A
              1  |               1       10       20  +         exon       A
              2  |               1       35       45  +         exon       A
              3  |               1       30       40  +         exon       A
              4  |               1        0       50  +         gene       B
              5  |               1       20       30  +         exon       B
        PyRanges with 6 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 1 strands.

        >>> gr.features.introns(feature_column="Feature", outer_feature="gene", inner_feature="exon", by="Id")
          index  |      Chromosome    Start      End  Strand    Feature    Id
          int64  |           int64    int64    int64  object    object     object
        -------  ---  ------------  -------  -------  --------  ---------  --------
              0  |               1        0       10  +         gene       A
              0  |               1       20       30  +         gene       A
              0  |               1       45      100  +         gene       A
              4  |               1        0       20  +         gene       B
              4  |               1       30       50  +         gene       B
        PyRanges with 5 rows, 6 columns, and 1 index columns (with 3 index duplicates).
        Contains 1 chromosomes and 1 strands.

        """
        gr = self.pr
        if gr.empty:
            return gr

        inner_df = mypy_ensure_pyranges(gr.loc[gr[feature_column] == inner_feature])
        outer_df = mypy_ensure_pyranges(gr.loc[gr[feature_column] == outer_feature])

        return outer_df.subtract_ranges(inner_df, match_by=by)


def _outside_bounds(df: DataFrame, **kwargs) -> DataFrame:
    df = df.copy()

    _chromsizes = kwargs.get("chromsizes")

    if isinstance(_chromsizes, pd.DataFrame):
        size_df = _chromsizes.df
        if not size_df.Chromosome.is_unique:
            msg = "Chromosomes must be unique in chromsizes."
            raise ValueError(msg)
        chromsizes = dict(*zip(size_df[CHROM_COL], size_df[END_COL], strict=True))
    else:
        if not isinstance(_chromsizes, dict):
            msg = "chromsizes must be a dictionary or a DataFrame"
            raise TypeError(msg)
        chromsizes = _chromsizes

    size = int(chromsizes[df.Chromosome.iloc[0]])
    clip = kwargs.get("clip", False)
    only_right = kwargs.get("only_right", False)

    ends_outright = df.End > size
    starts_outleft = df.Start < 0

    if not clip:  # i.e. remove
        df = df[~ends_outright] if only_right else df[~ends_outright & ~starts_outleft]

    else:
        starts_outright = df.Start >= size

        if only_right:
            df.loc[ends_outright, "End"] = size

            # removing intervals completely out of bounds
            df = df[~starts_outright]

        else:
            ends_outleft = df.End <= 0

            df.loc[ends_outright, "End"] = size
            df.loc[starts_outleft, "Start"] = 0

            # removing intervals completely out of bounds:
            df = df[~starts_outright & ~ends_outleft]

    return df


def genome_bounds(
    gr: "PyRanges",
    chromsizes: "dict[str | int, int] | PyRanges",
    *,
    clip: bool = False,
    only_right: bool = False,
) -> "PyRanges":
    """Remove or clip intervals outside of genome bounds.

    Parameters
    ----------
    gr : PyRanges
        Input intervals

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

    >>> pr.genomicfeatures.genome_bounds(gr, chromsizes)
      index  |      Chromosome    Start      End
      int64  |           int64    int64    int64
    -------  ---  ------------  -------  -------
          0  |               1        1        2
          1  |               3        5        7
    PyRanges with 2 rows, 3 columns, and 1 index columns.
    Contains 2 chromosomes.

    >>> pr.genomicfeatures.genome_bounds(gr, chromsizes, clip=True)
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

    >>> pr.genomicfeatures.genome_bounds(gr, chromsizes)
    Traceback (most recent call last):
    ...
    ValueError: Not all chromosomes were in the chromsize dict. This might mean that their types differed.
    Missing keys: {3}.
    Chromosome col had type: int64 while keys were of type: int

    """
    if isinstance(chromsizes, pd.DataFrame):
        chromsizes = dict(*zip(chromsizes[CHROM_COL], chromsizes[END_COL], strict=True))
    elif isinstance(chromsizes, dict):
        pass
    else:  # A hack because pyfaidx might not be installed, but we want type checking anyway
        pyfaidx_chromsizes = typing.cast(dict[str | int, list], chromsizes)
        chromsizes = {k: len(pyfaidx_chromsizes[k]) for k in pyfaidx_chromsizes}

    if missing_keys := set(gr[CHROM_COL]).difference(set(chromsizes.keys())):
        msg = f"""Not all chromosomes were in the chromsize dict. This might mean that their types differed.
Missing keys: {missing_keys}.
Chromosome col had type: {gr[CHROM_COL].dtype} while keys were of type: {', '.join({type(k).__name__ for k in chromsizes})}"""
        raise ValueError(msg)

    if not isinstance(chromsizes, dict):
        msg = "ERROR chromsizes must be a dictionary, or a PyRanges, or a pyfaidx.Fasta object"
        raise TypeError(msg)

    return (
        gr.groupby(CHROM_COL)
        .apply(
            _outside_bounds,
            chromsizes=chromsizes,
            clip=clip,
            only_right=only_right,
        )
        .reset_index(drop=True)
    )


def _last_tile(df: DataFrame, sizes: dict[str, int]) -> DataFrame:
    size = sizes[df.Chromosome.iloc[0]]
    df.iloc[-1, [*df.columns].index(END_COL)] = size
    return df


def tile_genome(
    chromsizes: "PyRanges | pd.DataFrame | dict[str | int, int]",
    tile_size: int,
    *,
    tile_last: bool = False,
) -> "PyRanges":
    """Create a tiled genome.

    Parameters
    ----------
    chromsizes : dict or PyRanges
        Dict or PyRanges describing the lengths of the chromosomes.

    tile_size : int
        Length of the tiles.

    tile_last : bool, default False
        Use chromosome length as end of last tile.

    See Also
    --------
    pyranges.PyRanges.tile : split intervals into adjacent non-overlapping tiles.

    Examples
    --------
    >>> import pyranges as pr
    >>> chromsizes = pr.example_data.chromsizes
    >>> chromsizes
    index    |    Chromosome    Start    End
    int64    |    category      int64    int64
    -------  ---  ------------  -------  ---------
    0        |    chr1          0        249250621
    1        |    chr2          0        243199373
    2        |    chr3          0        198022430
    3        |    chr4          0        191154276
    ...      |    ...           ...      ...
    21       |    chr19         0        59128983
    22       |    chr22         0        51304566
    23       |    chr21         0        48129895
    24       |    chrM          0        16571
    PyRanges with 25 rows, 3 columns, and 1 index columns.
    Contains 25 chromosomes.

    >>> pr.genomicfeatures.tile_genome(chromsizes, int(1e6))
    index    |    Chromosome    Start     End
    int64    |    category      int64     int64
    -------  ---  ------------  --------  --------
    0        |    chr1          0         1000000
    1        |    chr1          1000000   2000000
    2        |    chr1          2000000   3000000
    3        |    chr1          3000000   4000000
    ...      |    ...           ...       ...
    3110     |    chrY          56000000  57000000
    3111     |    chrY          57000000  58000000
    3112     |    chrY          58000000  59000000
    3113     |    chrY          59000000  59373566
    PyRanges with 3114 rows, 3 columns, and 1 index columns.
    Contains 25 chromosomes.

    """
    if isinstance(chromsizes, dict):
        chromsize_dict = chromsizes
        chromosomes, ends = list(chromsizes.keys()), list(chromsizes.values())
        df = pd.DataFrame({CHROM_COL: chromosomes, START_COL: 0, END_COL: ends})
        chromsizes = pd.DataFrame(df)
    else:
        chromsize_dict = dict(zip(chromsizes[CHROM_COL], chromsizes[END_COL], strict=True))

    gr = mypy_ensure_pyranges(chromsizes).tile(tile_size)

    if not tile_last:
        gr = gr.groupby(CHROM_COL).apply(_last_tile, sizes=chromsize_dict).reset_index(drop=True)

    return gr


def _keep_transcript_with_most_exons(df: pd.DataFrame) -> DataFrame:
    transcripts_with_most_exons = []

    for _, gdf in df.groupby("gene_id"):
        max_exon = gdf.exon_number.max()
        max_transcript = gdf.loc[gdf.exon_number == max_exon].Transcript.iloc[0]

        max_rows = gdf.loc[gdf.Transcript == max_transcript]

        transcripts_with_most_exons.append(max_rows)

    return pd.concat(transcripts_with_most_exons).reset_index(drop=True)


def filter_transcripts(df: pd.DataFrame) -> DataFrame:
    return _keep_transcript_with_most_exons(df)


def _tss(df: DataFrame, slack: int = 0) -> DataFrame:
    tss_pos = df.loc[df.Strand == "+"]

    tss_neg = df.loc[df.Strand == "-"].copy()

    tss_neg.loc[:, "Start"] = tss_neg.End - 1

    tss = pd.concat([tss_pos, tss_neg], sort=False)
    tss["End"] = tss.Start + 1
    tss.End = tss.End + slack
    tss.Start = tss.Start - slack
    tss.loc[tss.Start < 0, "Start"] = 0

    tss.index = pd.Index(range(len(tss)))

    return tss


def _tes(df: DataFrame, slack: int = 0) -> DataFrame:
    tes_pos = df.loc[df.Strand == "+"]

    tes_neg = df.loc[df.Strand == "-"].copy()

    tes_neg.loc[:, "End"] = tes_neg.Start + 1

    tes = pd.concat([tes_pos, tes_neg], sort=False)
    tes["Start"] = tes.End - 1
    tes.End = tes.End + slack
    tes.Start = tes.Start - slack
    tes.loc[tes.Start < 0, "Start"] = 0

    tes.index = pd.Index(range(len(tes)))

    return tes


by_to_id = {"gene": "gene_id", "transcript": "transcript_id"}
