from typing import Dict

import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame
from sorted_nearest.src.introns import find_introns  # type: ignore

import pyranges as pr
from pyranges.names import CHROM_COL, END_COL
from pyranges.pyranges_main import PyRanges

__all__ = ["genome_bounds", "tile_genome", "GenomicFeaturesMethods"]


class GenomicFeaturesMethods:

    """Namespace for methods using feature information.

    Accessed through `gr.features`."""

    def __init__(self, pr: PyRanges) -> None:
        self.pr = pr

    def tss(self) -> PyRanges:
        """Return the transcription start sites.

        Returns the 5' for every interval with feature "transcript".

        See Also
        --------
        pyranges.genomicfeatures.GenomicFeaturesMethods.tes : return the transcription end sites

        Examples
        --------

        >>> gr = pr.data.ensembl_gtf.get_with_loc_columns(["Source", "Feature"])
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
            raise Exception(
                "Cannot compute TSSes or TESes without strand info. Perhaps use extend() or subsequence() or spliced_subsequence() instead?"
            )

        gr = gr[gr.Feature == "transcript"]
        gr = gr.groupby(CHROM_COL).apply(_tss).reset_index(drop=True)

        gr.Feature = "tss"

        return pr.PyRanges(gr)

    def tes(self) -> PyRanges:
        """Return the transcription end sites.

        Returns the 3' for every interval with feature "transcript".

        See Also
        --------
        pyranges.genomicfeatures.GenomicFeaturesMethods.tss : return the transcription start sites

        Examples
        --------

        >>> gr = pr.data.ensembl_gtf.get_with_loc_columns(["Source", "Feature"])
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
            raise Exception(
                "Cannot compute TSSes or TESes without strand info. Perhaps use extend() or subsequence() or spliced_subsequence() instead?"
            )

        gr = gr[gr.Feature == "transcript"]
        gr = gr.groupby(CHROM_COL).apply(_tes).reset_index(drop=True)

        gr.Feature = "tes"

        return pr.PyRanges(gr)

    def introns(self, by: str = "gene") -> PyRanges:
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

        >>> gr = pr.data.ensembl_gtf.get_with_loc_columns(["Feature", "gene_id", "transcript_id"])
        >>> gr
        Chromosome    Start    End      Strand      Feature     gene_id          ...
        category      int64    int64    category    category    object           ...
        ------------  -------  -------  ----------  ----------  ---------------  -----
        1             11868    14409    +           gene        ENSG00000223972  ...
        1             11868    14409    +           transcript  ENSG00000223972  ...
        1             11868    12227    +           exon        ENSG00000223972  ...
        1             12612    12721    +           exon        ENSG00000223972  ...
        ...           ...      ...      ...         ...         ...              ...
        1             120724   133723   -           transcript  ENSG00000238009  ...
        1             133373   133723   -           exon        ENSG00000238009  ...
        1             129054   129223   -           exon        ENSG00000238009  ...
        1             120873   120932   -           exon        ENSG00000238009  ...
        PyRanges with 11 rows and 7 columns (1 columns not shown: "transcript_id").
        Contains 1 chromosomes and 2 strands.

        >>> gr.features.introns(by="gene")
        +--------------+------------+-----------+-----------+--------------+-----------------+-----------------+
        | Chromosome   | Feature    | Start     | End       | Strand       | gene_id         | transcript_id   |
        | (category)   | (object)   | (int64)   | (int64)   | (category)   | (object)        | (object)        |
        |--------------+------------+-----------+-----------+--------------+-----------------+-----------------|
        | 1            | intron     | 1173926   | 1174265   | +            | ENSG00000162571 | nan             |
        | 1            | intron     | 1174321   | 1174423   | +            | ENSG00000162571 | nan             |
        | 1            | intron     | 1174489   | 1174520   | +            | ENSG00000162571 | nan             |
        | 1            | intron     | 1175034   | 1179188   | +            | ENSG00000162571 | nan             |
        | ...          | ...        | ...       | ...       | ...          | ...             | ...             |
        | 1            | intron     | 874591    | 875046    | -            | ENSG00000283040 | nan             |
        | 1            | intron     | 875155    | 875525    | -            | ENSG00000283040 | nan             |
        | 1            | intron     | 875625    | 876526    | -            | ENSG00000283040 | nan             |
        | 1            | intron     | 876611    | 876754    | -            | ENSG00000283040 | nan             |
        +--------------+------------+-----------+-----------+--------------+-----------------+-----------------+
        Stranded PyRanges object has 311 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.

        >>> gr.features.introns(by="transcript")
        +--------------+------------+-----------+-----------+--------------+-----------------+-----------------+
        | Chromosome   | Feature    | Start     | End       | Strand       | gene_id         | transcript_id   |
        | (category)   | (object)   | (int64)   | (int64)   | (category)   | (object)        | (object)        |
        |--------------+------------+-----------+-----------+--------------+-----------------+-----------------|
        | 1            | intron     | 818202    | 818722    | +            | ENSG00000177757 | ENST00000326734 |
        | 1            | intron     | 960800    | 961292    | +            | ENSG00000187961 | ENST00000338591 |
        | 1            | intron     | 961552    | 961628    | +            | ENSG00000187961 | ENST00000338591 |
        | 1            | intron     | 961750    | 961825    | +            | ENSG00000187961 | ENST00000338591 |
        | ...          | ...        | ...       | ...       | ...          | ...             | ...             |
        | 1            | intron     | 732207    | 732980    | -            | ENSG00000230021 | ENST00000648019 |
        | 1            | intron     | 168165    | 169048    | -            | ENSG00000241860 | ENST00000655252 |
        | 1            | intron     | 165942    | 167958    | -            | ENSG00000241860 | ENST00000662089 |
        | 1            | intron     | 168165    | 169048    | -            | ENSG00000241860 | ENST00000662089 |
        +--------------+------------+-----------+-----------+--------------+-----------------+-----------------+
        Stranded PyRanges object has 1,043 rows and 7 columns from 1 chromosomes.
        For printing, the PyRanges was sorted on Chromosome and Strand.
        """

        assert by in ["gene", "transcript"]

        id_column = by_to_id[by]
        gr = self.pr.sort(id_column)

        if not len(gr):
            return pr.PyRanges()

        exons = gr.subset(lambda df: df.Feature == "exon")
        exons = exons.merge_overlaps(by=id_column)

        by_gr = gr.subset(lambda df: df.Feature == by)

        result = pyrange_apply(_introns2, by_gr, exons, **kwargs)

        return pr.from_dfs(result)


def _outside_bounds(df: DataFrame, **kwargs) -> DataFrame:
    df = df.copy()

    _chromsizes = kwargs.get("chromsizes")

    if isinstance(_chromsizes, PyRanges):
        size_df = _chromsizes.df
        if not size_df.Chromosome.is_unique:
            raise ValueError("Chromosomes must be unique in chromsizes.")
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
        if only_right:
            df = df[~ends_outright]
        else:
            df = df[~ends_outright & ~starts_outleft]

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
    gr: PyRanges,
    chromsizes: Dict[str, int],
    clip: bool = False,
    only_right: bool = False,
) -> PyRanges:
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

    if isinstance(chromsizes, pr.PyRanges):
        chromsizes = {k: v for k, v in zip(chromsizes.Chromosome, chromsizes.End)}

    elif isinstance(chromsizes, dict):
        pass

    else:
        try:
            import pyfaidx  # type: ignore

            if isinstance(chromsizes, pyfaidx.Fasta):
                chromsizes = {k: len(chromsizes[k]) for k in chromsizes.keys()}
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

    return pr.PyRanges(
        gr.groupby(CHROM_COL).apply(
            _outside_bounds, chromsizes=chromsizes, clip=clip, only_right=only_right,
        ).reset_index(drop=True)
    )


def _last_tile(df: DataFrame, sizes: dict[str, int]) -> DataFrame:
    size = sizes[df.Chromosome.iloc[0]]
    df.iloc[-1, df.columns.get_loc(END_COL)] = size
    return df


def tile_genome(
    chromsizes: PyRanges, tile_size: int, tile_last: bool = False,
) -> PyRanges:
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

    >>> chromsizes = pr.data.chromsizes
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
        df = pd.DataFrame({"Chromosome": chromosomes, "Start": 0, "End": ends})
        chromsizes = pr.PyRanges(df)
    else:
        chromsize_dict = dict(zip(chromsizes[CHROM_COL], chromsizes[END_COL], strict=True))

    gr = chromsizes.tile(tile_size)

    if not tile_last:
        gr = gr.groupby(CHROM_COL).apply(_last_tile, sizes=chromsize_dict).reset_index(drop=True)

    return pr.PyRanges(gr)


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

    # pd.options.mode.chained_assignment = None
    tss_neg.loc[:, "Start"] = tss_neg.End - 1

    # pd.options.mode.chained_assignment = "warn"
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

    # pd.options.mode.chained_assignment = None
    tes_neg.loc[:, "End"] = tes_neg.Start + 1

    # pd.options.mode.chained_assignment = "warn"
    tes = pd.concat([tes_pos, tes_neg], sort=False)
    tes["Start"] = tes.End - 1
    tes.End = tes.End + slack
    tes.Start = tes.Start - slack
    tes.loc[tes.Start < 0, "Start"] = 0

    tes.index = pd.Index(range(len(tes)))

    return tes


by_to_id = {"gene": "gene_id", "transcript": "transcript_id"}


def _introns2(df: DataFrame, exons: DataFrame, **kwargs) -> DataFrame:
    """TODO: refactor"""

    if df.empty or exons.empty:
        return pd.DataFrame(columns=df.columns)

    original_order = df.columns
    by = kwargs["by"]
    id_column = by_to_id[by]

    exons = exons[["Start", "End", id_column]]
    genes = df[["Start", "End", id_column]]
    exons.columns = pd.Index(["Start", "End", "by_id"])
    genes.columns = pd.Index(["Start", "End", "by_id"])

    intersection = pd.Series(np.intersect1d(exons["by_id"], genes["by_id"]))
    if len(intersection) == 0:
        return pd.DataFrame(columns=df.columns)

    exons = (
        exons[exons["by_id"].isin(intersection)]
        .reset_index(drop=True)
        .sort_values(["by_id", "Start"])
    )
    genes = (
        genes[genes["by_id"].isin(intersection)]
        .reset_index(drop=True)
        .sort_values(["by_id", "Start"])
    )
    df = df[df[id_column].isin(intersection)].reset_index(drop=True)

    assert len(genes) == len(
        genes.drop_duplicates("by_id")
    ), "The {id_column}s need to be unique to compute the introns.".format(
        id_column=id_column
    )

    exon_ids = exons["by_id"].shift() != exons["by_id"]
    by_ids = pd.Series(range(1, len(genes) + 1))
    df.insert(0, "__temp__", by_ids)

    if len(exons) > 1 and exons["by_id"].iloc[0] == exons["by_id"].iloc[1]:
        exon_ids.iloc[0] = False
        exon_ids = exon_ids.cumsum() + 1
    else:
        exon_ids = exon_ids.cumsum()

    assert (by_ids == exon_ids.drop_duplicates().values).all()
    starts, ends, ids = find_introns(
        genes.Start.values,
        genes.End.values,
        by_ids.values,
        exons.Start.values,
        exons.End.values,
        exon_ids.values,
    )

    introns = pd.DataFrame(
        data={
            "Chromosome": df.Chromosome.iloc[0],
            "Start": starts,
            "End": ends,
            "by_id": ids,
        }
    )

    vc = introns["by_id"].value_counts(sort=False).to_frame().reset_index()
    vc.columns = pd.Index(["by_id", "counts"])

    genes_without_introns = pd.DataFrame(
        data={
            "by_id": np.setdiff1d(np.array(by_ids.values), np.array(vc.by_id.values)),
            "counts": 0,
        }
    )

    vc = pd.concat([vc, genes_without_introns]).sort_values("by_id")

    original_ids = pd.Series(np.repeat(vc.by_id, vc.counts)).to_frame()
    original_ids = original_ids.merge(
        df[["__temp__", id_column]],
        right_on="__temp__",
        left_on="by_id",
        suffixes=("_drop", ""),
    )
    original_ids = original_ids.drop(
        ["__temp__"] + [c for c in original_ids.columns if c.endswith("_drop")], axis=1
    ).sort_values("by_id")
    introns.loc[:, "by_id"] = original_ids[id_column].values
    introns = introns.merge(
        df, left_on="by_id", right_on=id_column, suffixes=("", "_dropme")
    )
    introns = introns.drop(
        [c for c in introns.columns if c.endswith("_dropme")], axis=1
    )

    if (
        introns.Feature.dtype.name == "category"
        and "intron" not in introns.Feature.cat.categories
    ):
        introns.Feature.cat.add_categories(["intron"])
    introns.loc[:, "Feature"] = "intron"

    introns = introns[original_order]

    return introns
