import logging
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from natsort import natsorted  # type: ignore[import]

from pyranges.core.pyranges_helpers import ensure_pyranges

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from pyranges.core.pyranges_main import PyRanges

logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


def from_string(s: str) -> "PyRanges":
    """Create a PyRanges from multiline string.

    Parameters
    ----------
    s : str
        String with data.

    Examples
    --------
    >>> import pyranges as pr
    >>> s = '''Chromosome      Start        End Strand
    ... chr1  246719402  246719502      +
    ... chr5   15400908   15401008      +
    ... chr9   68366534   68366634      +
    ... chr14   79220091   79220191      +
    ... chr14  103456471  103456571      -'''

    >>> pr.from_string(s)
      index  |    Chromosome        Start        End  Strand
      int64  |    object            int64      int64  object
    -------  ---  ------------  ---------  ---------  --------
          0  |    chr1          246719402  246719502  +
          1  |    chr5           15400908   15401008  +
          2  |    chr9           68366534   68366634  +
          3  |    chr14          79220091   79220191  +
          4  |    chr14         103456471  103456571  -
    PyRanges with 5 rows, 4 columns, and 1 index columns.
    Contains 4 chromosomes and 2 strands.

    """
    from io import StringIO

    df = pd.read_csv(StringIO(s), sep=r"\s+", index_col=None)

    return ensure_pyranges(df)


def read_bed(f: Path, /, nrows: int | None = None) -> "PyRanges":
    """Return bed file as PyRanges.

    This is a reader for files that follow the bed format. They can have from
    3-12 columns which will be named like so:

    Chromosome Start End Name Score Strand ThickStart ThickEnd ItemRGB
    BlockCount BlockSizes BlockStarts

    Parameters
    ----------
    f : str
        Path to bed file

    nrows : Optional int, default None
        Number of rows to return.

    Notes
    -----
    If you just want to create a PyRanges from a tab-delimited bed-like file,
    use `pr.PyRanges(pandas.read_table(f))` instead.

    Returns
    -------
    PyRanges


    Examples
    --------
    >>> import pyranges as pr
    >>> path = pr.example_data.files["aorta.bed"]
    >>> pr.read_bed(path, nrows=5)
      index  |    Chromosome      Start      End  Name        Score  Strand
      int64  |    category        int64    int64  object      int64  category
    -------  ---  ------------  -------  -------  --------  -------  ----------
          0  |    chr1             9916    10115  H3K27me3        5  -
          1  |    chr1             9939    10138  H3K27me3        7  +
          2  |    chr1             9951    10150  H3K27me3        8  -
          3  |    chr1             9953    10152  H3K27me3        5  +
          4  |    chr1             9978    10177  H3K27me3        7  -
    PyRanges with 5 rows, 6 columns, and 1 index columns.
    Contains 1 chromosomes and 2 strands.

    """
    columns = [
        "Chromosome",
        "Start",
        "End",
        "Name",
        "Score",
        "Strand",
        "ThickStart",
        "ThickEnd",
        "ItemRGB",
        "BlockCount",
        "BlockSizes",
        "BlockStarts",
    ]
    path = Path(f)
    if path.name.endswith(".gz"):
        import gzip

        first_start = gzip.open(path).readline().decode().split()[1]  # noqa: SIM115
    else:
        first_start = path.open().readline().split()[1]

    header = None

    try:
        int(first_start)
    except ValueError:
        header = 0

    ncols = pd.read_table(path, nrows=2).shape[1]

    df = pd.read_csv(
        path,
        dtype={"Chromosome": "category", "Strand": "category"},
        nrows=nrows,
        header=header,
        names=columns[:ncols] if header != 0 else None,
        sep="\t",
    )

    df.columns = pd.Index(columns[: df.shape[1]])

    return ensure_pyranges(df)


def read_bam(
    f: str | Path,
    /,
    mapq: int = 0,
    required_flag: int = 0,
    filter_flag: int = 1540,
    *,
    sparse: bool = True,
) -> "PyRanges":
    """Return bam file as PyRanges.

    Parameters
    ----------
    f : str
        Path to bam file

    sparse : bool, default True
        Whether to return only the columns Chromosome, Start, End, Strand, Flag.
        Set to False to return also columns
        QueryStart, QueryEnd, QuerySequence, Name, Cigar, Quality (more time consuming).

    mapq : int, default 0
        Minimum mapping quality score.

    required_flag : int, default 0
        Flags which must be present for the interval to be read.

    filter_flag : int, default 1540
        Ignore reads with these flags. Default 1540, which means that either
        the read is unmapped, the read failed vendor or platfrom quality
        checks, or the read is a PCR or optical duplicate.

    Returns
    -------
    PyRanges


    Notes
    -----
    This functionality requires the library `bamread`. It can be installed with
    `pip install bamread` or `conda install -c bioconda bamread`.

    Examples
    --------
    >>> import pyranges as pr
    >>> path = pr.example_data.files["smaller.bam"]
    >>> pr.read_bam(path)
    index    |    Chromosome    Start     End       Strand      Flag
    int64    |    category      int64     int64     category    uint16
    -------  ---  ------------  --------  --------  ----------  --------
    0        |    chr1          887771    887796    -           16
    1        |    chr1          994660    994685    -           16
    2        |    chr1          1041102   1041127   +           0
    3        |    chr1          1770383   1770408   -           16
    ...      |    ...           ...       ...       ...         ...
    96       |    chr1          18800901  18800926  +           0
    97       |    chr1          18800901  18800926  +           0
    98       |    chr1          18855123  18855148  -           16
    99       |    chr1          19373470  19373495  +           0
    PyRanges with 100 rows, 5 columns, and 1 index columns.
    Contains 1 chromosomes and 2 strands.

    """
    path = Path(f)
    try:
        import bamread  # type: ignore[import]
    except ImportError:
        LOGGER.exception(
            "bamread must be installed to read bam. Use `conda install -c bioconda bamread` or `pip install bamread` to install it.",
        )
        sys.exit(1)

    if bamread.__version__ in {
        "0.0.1",
        "0.0.2",
        "0.0.3",
        "0.0.4",
        "0.0.5",
        "0.0.6",
        "0.0.7",
        "0.0.8",
        "0.0.9",
    }:
        LOGGER.exception(
            "bamread not recent enough. Must be 0.0.10 or higher. Use `conda install -c bioconda 'bamread>=0.0.10'` or `pip install bamread>=0.0.10` to install it.",
        )
        sys.exit(1)

    if sparse:
        return ensure_pyranges(bamread.read_bam(path, mapq, required_flag, filter_flag))
    df = bamread.read_bam_full(path, mapq, required_flag, filter_flag)
    return ensure_pyranges(df)


def _fetch_gene_transcript_exon_id(attribute: pd.Series, annotation: str | None = None) -> pd.DataFrame:
    no_quotes = attribute.str.replace('"', "").str.replace("'", "")

    df = no_quotes.str.extract(
        "gene_id.?(.+?);(?:.*group_by.?(.+?);)?(?:.*exon_number.?(.+?);)?(?:.*exon_id.?(.+?);)?",
        expand=True,
    )  # .iloc[:, [1, 2, 3]]

    df.columns = pd.Index(["gene_id", "group_by", "exon_number", "exon_id"])

    if annotation == "ensembl":
        newdfs = []
        for c in ["gene_id", "group_by", "exon_id"]:
            r = df[c].astype(str).str.extract(r"(\d+)").astype(float)
            newdfs.append(r)

        newdf = pd.concat(newdfs, axis=1)
        newdf.insert(2, "exon_number", df["exon_number"])
        df = newdf

    return df


def find_first_data_line_index(file_path: Path) -> int:
    """Find the first line that is not a comment."""
    first_non_comment = 0
    try:
        import gzip

        with gzip.open(file_path) as zh:
            for i, zl in enumerate(zh):
                if zl.decode()[0] != "#":
                    first_non_comment = i
                    break
    except (OSError, TypeError):  # not a gzipped file, or StringIO
        fh = file_path.open()
        for i, line in enumerate(fh):
            if line[0] != "#":
                first_non_comment = i
                break
        fh.close()

    return first_non_comment


def read_gtf(
    f: str | Path,
    /,
    *,
    nrows: bool | None = None,
    full: bool = True,
    duplicate_attr: bool = False,
    ignore_bad: bool = False,
) -> "PyRanges":
    r"""Read files in the Gene Transfer Format.

    Parameters
    ----------
    f : str
        Path to GTF file.

    full : bool, default True
        Whether to read and interpret the annotation column.

    nrows : int, default None
        Number of rows to read. Default None, i.e. all.

    duplicate_attr : bool, default False
        Whether to handle (potential) duplicate attributes or just keep last one.

    ignore_bad : bool, default False
        Whether to ignore bad lines or raise an error.


    Returns
    -------
    PyRanges

    Note
    ----
    The GTF format encodes both Start and End as 1-based included.
    PyRanges encodes intervals as 0-based, Start included and End excluded.

    See Also
    --------
    pyranges.read_gff3 : read files in the General Feature Format

    Examples
    --------
    >>> import pyranges as pr
    >>> from tempfile import NamedTemporaryFile
    >>> contents = ['#!genome-build GRCh38.p10']
    >>> contents.append('1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";')
    >>> contents.append('1\thavana\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; group_by "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; tag "basic"; transcript_support_level "1";')
    >>> f = NamedTemporaryFile("w")
    >>> _bytes_written = f.write("\n".join(contents))
    >>> f.flush()
    >>> pr.read_gtf(f.name)
      index  |      Chromosome  Source    Feature       Start      End  Score     Strand      Frame     gene_id          ...
      int64  |        category  object    category      int64    int64  object    category    object    object           ...
    -------  ---  ------------  --------  ----------  -------  -------  --------  ----------  --------  ---------------  -----
          0  |               1  havana    gene          11868    14409  .         +           .         ENSG00000223972  ...
          1  |               1  havana    transcript    11868    14409  .         +           .         ENSG00000223972  ...
    PyRanges with 2 rows, 20 columns, and 1 index columns. (11 columns not shown: "gene_version", "gene_name", "gene_source", ...).
    Contains 1 chromosomes and 1 strands.

    """
    path = Path(f)
    _skiprows = find_first_data_line_index(path)

    if full:
        gr = read_gtf_full(
            path,
            nrows=nrows,
            skiprows=_skiprows,
            duplicate_attr=duplicate_attr,
            ignore_bad=ignore_bad,
        )
    else:
        gr = read_gtf_restricted(path, _skiprows, nrows=None)

    return gr


def read_gtf_full(
    f: str | Path,
    /,
    nrows: int | None = None,
    skiprows: int = 0,
    chunksize: int = int(1e5),  # for unit-testing purposes
    *,
    duplicate_attr: bool = False,
    ignore_bad: bool = False,
) -> "PyRanges":
    """Read files in the Gene Transfer Format into a PyRanges, including the annotation column.

    Parameters
    ----------
    f : str
        Path to GTF file.

    nrows : int, default None
        Number of rows to read. Default None, i.e. all.

    skiprows : int, default 0
        Number of rows to skip. Default 0.

    chunksize : int, default 100000
        Number of rows to read at a time. Default 100000.

    duplicate_attr : bool, default False
        Whether to handle (potential) duplicate attributes or just keep last one.

    ignore_bad : bool, default False
        Whether to ignore bad lines or raise an error.

    Returns
    -------
    PyRanges

    """
    dtypes: Mapping = {"Chromosome": "category", "Feature": "category", "Strand": "category"}

    names = ["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute"]
    path = Path(f)

    df_iter = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=names,
        dtype=dtypes,
        chunksize=chunksize,
        skiprows=skiprows,
        nrows=nrows,
    )

    _to_rows = to_rows_keep_duplicates if duplicate_attr else to_rows

    dfs = []
    for df in df_iter:
        extra = _to_rows(df.Attribute.astype(str), ignore_bad=ignore_bad)
        _df = df.drop("Attribute", axis=1)
        extra = extra.set_index(_df.index)
        ndf = pd.concat([_df, extra], axis=1, sort=False)
        dfs.append(ndf)

    df = pd.concat(dfs, sort=False)
    df.loc[:, "Start"] = df.Start - 1

    return ensure_pyranges(df)


def parse_kv_fields(line: str) -> list[list[str]]:
    """Parse GTF attribute column."""
    return [kv.replace('""', '"NA"').replace('"', "").split(None, 1) for kv in line.rstrip("; ").split("; ")]


def to_rows(anno: pd.Series, *, ignore_bad: bool = False) -> pd.DataFrame:
    """Parse GTF attribute column into a dataframe of attribute columns."""
    entry = ""
    try:
        row = anno.head(1)
        for entry in row:
            str(entry).replace('"', "").replace(";", "").split()
    except AttributeError:
        msg = f"Invalid attribute string: {entry}. If the file is in GFF3 format, use pr.read_gff3 instead."
        raise AttributeError(msg) from AttributeError

    rowdicts = []
    line = ""
    try:
        for line in anno:
            rowdicts.append(dict(parse_kv_fields(line)))
    except ValueError:
        if not ignore_bad:
            LOGGER.exception(
                "The following line is not parseable as gtf:\n%s\n\nTo ignore bad lines use ignore_bad=True.",
                line,
            )
            raise

    return pd.DataFrame.from_records(rowdicts)


def to_rows_keep_duplicates(anno: pd.Series, *, ignore_bad: bool = False) -> pd.DataFrame:
    """If an entry is found multiple times in the attribute string, keep all of them.

    Examples
    --------
    >>> anno = pd.Series(["gene DDX11L1; gene sonic; unique hi"])
    >>> result = to_rows_keep_duplicates(anno)
    >>> result.to_dict(orient="records")
    [{'gene': 'DDX11L1,sonic', 'unique': 'hi'}]

    """
    rowdicts = []
    line = ""
    try:
        for line in anno:
            rowdict = {}

            # rstrip: allows for GFF not having a last ";", or having final spaces
            for k, v in tuple(parse_kv_fields(line)):
                if k not in rowdict:
                    rowdict[k] = [v]
                else:
                    rowdict[k].append(v)

            rowdicts.append({k: ",".join(v) if isinstance(v, list) else v for k, v in rowdict.items()})
    except ValueError:
        if not ignore_bad:
            LOGGER.exception(
                "The following line is not parseable as gtf:\n\n%s\n\nTo ignore bad lines use ignore_bad=True.",
                line,
            )

    return pd.DataFrame.from_records(rowdicts)


def read_gtf_restricted(f: str | Path, skiprows: int | None, nrows: int | None = None) -> "PyRanges":
    """Read certain columns from GTF file.

    Seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    # source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position of the feature, with sequence numbering starting at 1.
    end - End position of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    # frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

    """
    dtypes: Mapping = {"Chromosome": "category", "Feature": "category", "Strand": "category"}
    path = Path(f)
    usecols: Sequence[int] = [0, 2, 3, 4, 5, 6, 8]

    df_iter = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        usecols=np.array(usecols),
        header=None,
        names=["Chromosome", "Feature", "Start", "End", "Score", "Strand", "Attribute"],
        dtype=dtypes,
        chunksize=int(1e5),
        skiprows=skiprows if skiprows is not None else False,
        nrows=nrows,
    )

    dfs = []
    for df in df_iter:
        if sum(df.Score == ".") == len(df):
            cols_to_concat = ["Chromosome", "Start", "End", "Strand", "Feature"]
        else:
            cols_to_concat = ["Chromosome", "Start", "End", "Strand", "Feature", "Score"]

        extract = _fetch_gene_transcript_exon_id(df.Attribute)
        extract.columns = pd.Index(["gene_id", "group_by", "exon_number", "exon_id"])

        extract.exon_number = extract.exon_number.astype(float)

        extract = extract.set_index(df.index)
        _df = pd.concat([df[cols_to_concat], extract], axis=1, sort=False)

        dfs.append(_df)

    df = pd.concat(dfs, sort=False)

    df.loc[:, "Start"] = df.Start - 1

    return ensure_pyranges(df)


def to_rows_gff3(anno: pd.Series) -> pd.DataFrame:
    """Parse GFF3 attribute column into a dataframe of attribute columns."""
    rowdicts = [to_keys_and_values(line) for line in list(anno)]

    return pd.DataFrame.from_records(rowdicts).set_index(anno.index)


def to_keys_and_values(line: str) -> dict[str, str]:
    """Parse GFF3 attribute column."""
    return dict(it.split("=") for it in line.rstrip("; ").split(";"))


def read_gff3(
    f: str | Path,
    nrows: int | None = None,
    *,
    full: bool = True,
) -> "PyRanges":
    """Read files in the General Feature Format into a PyRanges.

    Parameters
    ----------
    f : str
        Path to GFF file.

    full : bool, default True
        Whether to read and interpret the annotation column.

    nrows : int, default None
        Number of rows to read. Default None, i.e. all.

    Returns
    -------
    PyRanges

    Notes
    -----
    The gff3 format encodes both Start and End as 1-based included.
    PyRanges (and also the DF returned by this function, if as_df=True), instead
    encodes intervals as 0-based, Start included and End excluded.

    See Also
    --------
    pyranges.read_gtf : read files in the Gene Transfer Format

    """
    path = Path(f)
    _skiprows = find_first_data_line_index(path)

    if not full:
        return read_gtf_restricted(path, _skiprows, nrows=nrows)

    dtypes: Mapping = {"Chromosome": "category", "Feature": "category", "Strand": "category"}

    names = ["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute"]

    df_iter = pd.read_csv(
        path,
        comment="#",
        sep="\t",
        header=None,
        names=names,
        dtype=dtypes,
        chunksize=int(1e5),
        skiprows=_skiprows,
        nrows=nrows,
    )

    dfs = []
    for df in df_iter:
        extra = to_rows_gff3(df.Attribute.astype(str))
        _df = df.drop("Attribute", axis=1)
        extra = extra.set_index(_df.index)
        ndf = pd.concat([_df, extra], axis=1, sort=False)
        dfs.append(ndf)

    df = pd.concat(dfs, sort=False)

    df.loc[:, "Start"] = df.Start - 1

    return ensure_pyranges(df)


def read_bigwig(f: str | Path) -> "PyRanges":
    """Read bigwig files into a PyRanges.

    Parameters
    ----------
    f : str
        Path to bw file.

    Returns
    -------
    PyRanges

    Note
    ----
    This function requires the library pyBigWig, it can be installed with pip install pyBigWig

    Examples
    --------
    >>> import pyranges as pr
    >>> path = pr.example_data.files["bigwig.bw"]
    >>> pr.read_bigwig(path)
      index  |      Chromosome    Start      End      Value
      int64  |          object    int64    int64    float64
    -------  ---  ------------  -------  -------  ---------
          0  |               1        0        1        0.1
          1  |               1        1        2        0.2
          2  |               1        2        3        0.3
          3  |               1      100      150        1.4
          4  |               1      150      151        1.5
          5  |              10      200      300        2
    PyRanges with 6 rows, 4 columns, and 1 index columns.
    Contains 2 chromosomes.


    """
    try:
        import pyBigWig  # type: ignore[import]
    except ModuleNotFoundError:
        LOGGER.exception(
            "pyBigWig must be installed to read bigwigs. Use `pip install pyBigWig` to install it.",
        )
        sys.exit(1)

    path = Path(f)
    bw = pyBigWig.open(str(path))

    size = int(1e5)
    chromosomes = bw.chroms()

    dfs = {}

    for chromosome in natsorted(chromosomes):
        outstarts = []
        outends = []
        outvalues = []

        length = chromosomes[chromosome]

        starts = list(range(0, length, size))
        ends = list(range(size, length + size, size))
        ends[-1] = length
        for start, end in zip(starts, ends, strict=True):
            intervals = bw.intervals(chromosome, start, end)
            if intervals is not None:
                for s, e, v in intervals:
                    outstarts.append(s)
                    outends.append(e)
                    outvalues.append(v)

        outstarts = pd.Series(outstarts)
        outends = pd.Series(outends)
        outvalues = pd.Series(outvalues)
        dfs[chromosome] = pd.DataFrame(
            {
                "Chromosome": chromosome,
                "Start": outstarts,
                "End": outends,
                "Value": outvalues,
            },
        )

    return ensure_pyranges(pd.concat(dfs).reset_index(drop=True))
