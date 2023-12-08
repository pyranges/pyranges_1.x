from typing import TYPE_CHECKING

import pandas as pd
from pandas.core.frame import DataFrame

from pyranges.names import CHROM_COL, END_COL, START_COL

__all__ = ["genome_bounds", "tile_genome", "GenomicFeaturesMethods"]

if TYPE_CHECKING:
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
        Chromosome    Start    End      Strand      Source    Feature
        category      int64    int64    category    object    category
        ------------  -------  -------  ----------  --------  ----------
        1             11868    14409    +           havana    gene
        1             11868    14409    +           havana    transcript
        1             11868    12227    +           havana    exon
        1             12612    12721    +           havana    exon
        ...           ...      ...      ...         ...       ...
        1             120724   133723   -           ensembl   transcript
        1             133373   133723   -           ensembl   exon
        1             129054   129223   -           ensembl   exon
        1             120873   120932   -           ensembl   exon
        PyRanges with 11 rows and 6 columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.features.tss()
          Chromosome    Start      End  Strand      Source    Feature
            category    int64    int64  category    object    object
        ------------  -------  -------  ----------  --------  ---------
                   1    11868    11869  +           havana    tss
                   1   133722   133723  -           ensembl   tss
        PyRanges with 2 rows and 6 columns.
        Contains 1 chromosomes and 2 strands.
        """
        gr = self.pr

        if not gr.strand_values_valid:
            msg = (
                "Cannot compute TSSes or TESes without strand info. Perhaps use extend()"
                "or subsequence() or spliced_subsequence() instead?"
            )
            raise AssertionError(msg)

        gr = gr[gr.Feature == "transcript"]
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
        Chromosome    Start    End      Strand      Source    Feature
        category      int64    int64    category    object    category
        ------------  -------  -------  ----------  --------  ----------
        1             11868    14409    +           havana    gene
        1             11868    14409    +           havana    transcript
        1             11868    12227    +           havana    exon
        1             12612    12721    +           havana    exon
        ...           ...      ...      ...         ...       ...
        1             120724   133723   -           ensembl   transcript
        1             133373   133723   -           ensembl   exon
        1             129054   129223   -           ensembl   exon
        1             120873   120932   -           ensembl   exon
        PyRanges with 11 rows and 6 columns.
        Contains 1 chromosomes and 2 strands.

        >>> gr.features.tes()
          Chromosome    Start      End  Strand      Source    Feature
            category    int64    int64  category    object    object
        ------------  -------  -------  ----------  --------  ---------
                   1    14408    14409  +           havana    tes
                   1   120724   120725  -           ensembl   tes
        PyRanges with 2 rows and 6 columns.
        Contains 1 chromosomes and 2 strands.
        """
        gr = self.pr

        if not gr.strand_values_valid:
            msg = "Cannot compute TSSes or TESes without strand info. Perhaps use extend() or subsequence() or spliced_subsequence() instead?"
            raise ValueError(msg)

        gr = gr[gr.Feature == "transcript"]
        gr = gr.groupby(CHROM_COL).apply(_tes).reset_index(drop=True)

        gr.Feature = "tes"

        return gr

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
          Chromosome    Start      End  Strand      Feature     gene_id          transcript_id
            category    int64    int64  category    category    object           object
        ------------  -------  -------  ----------  ----------  ---------------  ---------------
                   1    11868    14409  +           gene        ENSG00000223972  nan
                   1    11868    14409  +           transcript  ENSG00000223972  ENST00000456328
                   1    11868    12227  +           exon        ENSG00000223972  ENST00000456328
                   1    12612    12721  +           exon        ENSG00000223972  ENST00000456328
                   1    13220    14409  +           exon        ENSG00000223972  ENST00000456328
        PyRanges with 5 rows and 7 columns.
        Contains 1 chromosomes and 1 strands.

        >>> gr = pr.from_string('''Chromosome   Start      End  Strand      Feature Id
        ...          1   0    100  +           gene A
        ...          1   10    20  +           exon A
        ...          1   35    45  +           exon A
        ...          1   30    40  +           exon A
        ...          1   0    50   +           gene B
        ...          1   20   30  +           exon  B''')
        >>> gr
          Chromosome    Start      End  Strand    Feature    Id
               int64    int64    int64  object    object     object
        ------------  -------  -------  --------  ---------  --------
                   1        0      100  +         gene       A
                   1       10       20  +         exon       A
                   1       35       45  +         exon       A
                   1       30       40  +         exon       A
                   1        0       50  +         gene       B
                   1       20       30  +         exon       B
        PyRanges with 6 rows and 6 columns.
        Contains 1 chromosomes and 1 strands.
        >>> gr.features.introns(feature_column="Feature", outer_feature="gene", inner_feature="exon", by="Id")
          Chromosome    Start      End  Strand    Feature    Id
               int64    int64    int64  object    object     object
        ------------  -------  -------  --------  ---------  --------
                   1        0       10  +         gene       A
                   1       20       30  +         gene       A
                   1       45      100  +         gene       A
                   1        0       20  +         gene       B
                   1       30       50  +         gene       B
        PyRanges with 5 rows and 6 columns.
        Contains 1 chromosomes and 1 strands.
        """
        gr = self.pr
        if gr.empty:
            return gr

        inner_df = gr[gr[feature_column] == inner_feature]
        outer_df = gr[gr[feature_column] == outer_feature]

        return outer_df.subtract_intervals(inner_df, by=by)


def _outside_bounds(df: DataFrame, **kwargs) -> DataFrame:
    df = df.copy()

    _chromsizes = kwargs.get("chromsizes")

    if isinstance(_chromsizes, pd.DataFrame):
        size_df = _chromsizes.df
        if not size_df.Chromosome.is_unique:
            msg = "Chromosomes must be unique in chromsizes."
            raise ValueError(msg)
        chromsizes = {k: v for k, v in zip(size_df.Chromosome, size_df.End)}
    else:
        assert isinstance(_chromsizes, dict)
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
    chromsizes: dict[str, int],
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
      Chromosome      Start        End
           int64      int64      int64
    ------------  ---------  ---------
               1          1          2
               1  249250600  249250640
               3          5          7
    PyRanges with 3 rows and 3 columns.
    Contains 2 chromosomes.

    >>> chromsizes = {1: 249250621, 3: 500}
    >>> chromsizes
    {1: 249250621, 3: 500}

    >>> pr.gf.genome_bounds(gr, chromsizes)
      Chromosome    Start      End
           int64    int64    int64
    ------------  -------  -------
               1        1        2
               3        5        7
    PyRanges with 2 rows and 3 columns.
    Contains 2 chromosomes.

    >>> pr.gf.genome_bounds(gr, chromsizes, clip=True)
      Chromosome      Start        End
           int64      int64      int64
    ------------  ---------  ---------
               1          1          2
               1  249250600  249250621
               3          5          7
    PyRanges with 3 rows and 3 columns.
    Contains 2 chromosomes.

    >>> del chromsizes[3]
    >>> chromsizes
    {1: 249250621}

    >>> pr.gf.genome_bounds(gr, chromsizes)
    Traceback (most recent call last):
    ...
    ValueError: Not all chromosomes were in the chromsize dict. This might mean that their types differed.
    Missing keys: {3}.
    Chromosome col had type: int64 while keys were of type: int
    """
    if isinstance(chromsizes, pd.DataFrame):
        chromsizes = {k: v for k, v in zip(chromsizes.Chromosome, chromsizes.End)}

    elif isinstance(chromsizes, dict):
        pass

    else:
        try:
            import pyfaidx  # type: ignore[import]

            if isinstance(chromsizes, pyfaidx.Fasta):
                chromsizes = {k: len(chromsizes[k]) for k in chromsizes}
        except ImportError:
            pass

    if missing_keys := set(gr[CHROM_COL]).difference(set(chromsizes.keys())):
        msg = f"""Not all chromosomes were in the chromsize dict. This might mean that their types differed.
Missing keys: {missing_keys}.
Chromosome col had type: {gr[CHROM_COL].dtype} while keys were of type: {', '.join(set(type(k).__name__ for k in chromsizes))}"""
        raise ValueError(msg)

    assert isinstance(
        chromsizes, dict
    ), "ERROR chromsizes must be a dictionary, or a PyRanges, or a pyfaidx.Fasta object"

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
    df.iloc[-1, df.columns.get_loc(END_COL)] = size
    return df


def tile_genome(
    chromsizes: "PyRanges",
    tile_size: int,
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
    Chromosome    Start    End
    category      int64    int64
    ------------  -------  ---------
    chr1          0        249250621
    chr2          0        243199373
    chr3          0        198022430
    chr4          0        191154276
    ...           ...      ...
    chr19         0        59128983
    chr22         0        51304566
    chr21         0        48129895
    chrM          0        16571
    PyRanges with 25 rows and 3 columns.
    Contains 25 chromosomes.

    >>> pr.gf.tile_genome(chromsizes, int(1e6))
    Chromosome    Start     End
    category      int64     int64
    ------------  --------  --------
    chr1          0         1000000
    chr1          1000000   2000000
    chr1          2000000   3000000
    chr1          3000000   4000000
    ...           ...       ...
    chrY          56000000  57000000
    chrY          57000000  58000000
    chrY          58000000  59000000
    chrY          59000000  59373566
    PyRanges with 3114 rows and 3 columns.
    Contains 25 chromosomes.
    """
    if isinstance(chromsizes, dict):
        chromsize_dict = chromsizes
        chromosomes, ends = list(chromsizes.keys()), list(chromsizes.values())
        df = pd.DataFrame({CHROM_COL: chromosomes, START_COL: 0, END_COL: ends})
        chromsizes = pd.DataFrame(df)
    else:
        chromsize_dict = dict(zip(chromsizes[CHROM_COL], chromsizes[END_COL], strict=True))

    gr = chromsizes.tile(tile_size)

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
