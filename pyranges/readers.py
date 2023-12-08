import sys
from pathlib import Path

import pandas as pd
from natsort import natsorted  # type: ignore[import]

from pyranges.pyranges_main import PyRanges


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
    Chromosome        Start        End  Strand
    object            int64      int64  object
    ------------  ---------  ---------  --------
    chr1          246719402  246719502  +
    chr5           15400908   15401008  +
    chr9           68366534   68366634  +
    chr14          79220091   79220191  +
    chr14         103456471  103456571  -
    PyRanges with 5 rows and 4 columns.
    Contains 4 chromosomes and 2 strands.
    """
    from io import StringIO

    df = pd.read_csv(StringIO(s), sep=r"\s+", index_col=None)

    return PyRanges(df)


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

    Examples
    --------
    >>> import pyranges as pr
    >>> path = pr.example_data.files["aorta.bed"]
    >>> pr.read_bed(path, nrows=5)
    Chromosome      Start      End  Name        Score  Strand
    category        int64    int64  object      int64  category
    ------------  -------  -------  --------  -------  ----------
    chr1             9916    10115  H3K27me3        5  -
    chr1             9939    10138  H3K27me3        7  +
    chr1             9951    10150  H3K27me3        8  -
    chr1             9953    10152  H3K27me3        5  +
    chr1             9978    10177  H3K27me3        7  -
    PyRanges with 5 rows and 6 columns.
    Contains 1 chromosomes and 2 strands.
    """
    columns = (
        "Chromosome Start End Name Score Strand ThickStart ThickEnd ItemRGB BlockCount BlockSizes BlockStarts".split()
    )
    path = Path(f)
    if path.name.endswith(".gz"):
        import gzip

        first_start = gzip.open(path).readline().decode().split()[1]
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

    return PyRanges(df)


def read_bam(
    f: str | Path,
    /,
    sparse: bool = True,
    mapq: int = 0,
    required_flag: int = 0,
    filter_flag: int = 1540,
) -> "PyRanges":
    """Return bam file as PyRanges.

    Parameters
    ----------
    f : str

        Path to bam file

    sparse : bool, default True

        Whether to return only.

    mapq : int, default 0

        Minimum mapping quality score.

    required_flag : int, default 0

        Flags which must be present for the interval to be read.

    filter_flag : int, default 1540

        Ignore reads with these flags. Default 1540, which means that either
        the read is unmapped, the read failed vendor or platfrom quality
        checks, or the read is a PCR or optical duplicate.

    Notes
    -----
    This functionality requires the library `bamread`. It can be installed with
    `pip install bamread` or `conda install -c bioconda bamread`.

    Examples
    --------
    >>> import pyranges as pr
    >>> path = pr.example_data.files["smaller.bam"]
    >>> pr.read_bam(path)
    Chromosome    Start     End       Strand      Flag
    category      int64     int64     category    uint16
    ------------  --------  --------  ----------  --------
    chr1          887771    887796    -           16
    chr1          994660    994685    -           16
    chr1          1041102   1041127   +           0
    chr1          1770383   1770408   -           16
    ...           ...       ...       ...         ...
    chr1          18800901  18800926  +           0
    chr1          18800901  18800926  +           0
    chr1          18855123  18855148  -           16
    chr1          19373470  19373495  +           0
    PyRanges with 100 rows and 5 columns.
    Contains 1 chromosomes and 2 strands.
    """
    path = Path(f)
    try:
        import bamread  # type: ignore[import]
    except ImportError:
        print(
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
        print(
            "bamread not recent enough. Must be 0.0.10 or higher. Use `conda install -c bioconda 'bamread>=0.0.10'` or `pip install bamread>=0.0.10` to install it.",
        )
        sys.exit(1)

    if sparse:
        df = bamread.read_bam(path, mapq, required_flag, filter_flag)
    else:
        try:
            df = bamread.read_bam_full(path, mapq, required_flag, filter_flag)
        except AttributeError:
            print("bamread version 0.0.6 or higher is required to read bam non-sparsely.")

    return PyRanges(df)


def _fetch_gene_transcript_exon_id(attribute: pd.Series, annotation: str | None = None) -> pd.DataFrame:
    no_quotes = attribute.str.replace('"', "").str.replace("'", "")

    df = no_quotes.str.extract(
        "gene_id.?(.+?);(?:.*transcript_id.?(.+?);)?(?:.*exon_number.?(.+?);)?(?:.*exon_id.?(.+?);)?",
        expand=True,
    )  # .iloc[:, [1, 2, 3]]

    df.columns = pd.Index("gene_id transcript_id exon_number exon_id".split())

    if annotation == "ensembl":
        newdfs = []
        for c in "gene_id transcript_id exon_id".split():
            r = df[c].astype(str).str.extract(r"(\d+)").astype(float)
            newdfs.append(r)

        newdf = pd.concat(newdfs, axis=1)
        newdf.insert(2, "exon_number", df["exon_number"])
        df = newdf

    return df


def skiprows(f: Path) -> int:
    try:
        import gzip

        zh = gzip.open(f)
        for i, zl in enumerate(zh):
            if zl.decode()[0] != "#":
                break
        zh.close()
    except (OSError, TypeError):  # not a gzipped file, or StringIO
        fh = f.open()
        for i, line in enumerate(fh):
            if line[0] != "#":
                break
        fh.close()

    return i


def read_gtf(
    f: str | Path,
    /,
    full: bool = True,
    nrows: bool | None = None,
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
    >>> contents.append('1\thavana\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; tag "basic"; transcript_support_level "1";')
    >>> f = NamedTemporaryFile("w")
    >>> _bytes_written = f.write("\n".join(contents))
    >>> f.flush()
    >>> pr.read_gtf(f.name)
      Chromosome  Source    Feature       Start      End  Score     Strand      Frame     gene_id          ...
        category  object    category      int64    int64  object    category    object    object           ...
    ------------  --------  ----------  -------  -------  --------  ----------  --------  ---------------  -----
               1  havana    gene          11868    14409  .         +           .         ENSG00000223972  ...
               1  havana    transcript    11868    14409  .         +           .         ENSG00000223972  ...
    PyRanges with 2 rows and 20 columns (11 columns not shown: "gene_version", "gene_name", "gene_source", ...).
    Contains 1 chromosomes and 1 strands.
    """
    path = Path(f)
    _skiprows = skiprows(path)

    if full:
        gr = read_gtf_full(
            path,
            nrows,
            _skiprows,
            duplicate_attr,
            ignore_bad=ignore_bad,
        )
    else:
        gr = read_gtf_restricted(path, _skiprows, nrows=None)

    return gr


def read_gtf_full(
    f: str | Path,
    nrows: int | None = None,
    skiprows: int = 0,
    duplicate_attr: bool = False,
    ignore_bad: bool = False,
    chunksize: int = int(1e5),  # for unit-testing purposes
) -> "PyRanges":
    dtypes = {"Chromosome": "category", "Feature": "category", "Strand": "category"}

    names = "Chromosome Source Feature Start End Score Strand Frame Attribute".split()
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

    return PyRanges(df)


def parse_kv_fields(line: str) -> list[list[str]]:
    # rstrip: allows for GFF not having a last ";", or having final spaces
    return [kv.replace('""', '"NA"').replace('"', "").split(None, 1) for kv in line.rstrip("; ").split("; ")]


def to_rows(anno: pd.Series, ignore_bad: bool = False) -> pd.DataFrame:
    try:
        row = anno.head(1)
        for entry in row:
            str(entry).replace('"', "").replace(";", "").split()
    except AttributeError:
        msg = f"Invalid attribute string: {entry}. If the file is in GFF3 format, use pr.read_gff3 instead."
        raise AttributeError(msg) from AttributeError

    rowdicts = []
    try:
        for line in anno:
            rowdicts.append({k: v for k, v in parse_kv_fields(line)})  # noqa: PERF401
    except ValueError:
        if not ignore_bad:
            print(f"The following line is not parseable as gtf:\n{line}\n\nTo ignore bad lines use ignore_bad=True.")
            raise

    return pd.DataFrame.from_records(rowdicts)


def to_rows_keep_duplicates(anno: pd.Series, ignore_bad: bool = False) -> pd.DataFrame:
    """If an entry is found multiple times in the attribute string, keep all of them.

    Examples
    --------
    >>> anno = pd.Series(["gene DDX11L1; gene sonic; unique hi"])
    >>> result = to_rows_keep_duplicates(anno)
    >>> result.to_dict(orient="records")
    [{'gene': 'DDX11L1,sonic', 'unique': 'hi'}]
    """
    rowdicts = []
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
            print(f"The following line is not parseable as gtf:\n\n{line}\n\nTo ignore bad lines use ignore_bad=True.")
            raise

    return pd.DataFrame.from_records(rowdicts)


def read_gtf_restricted(f: str | Path, skiprows: int | None, nrows: int | None = None) -> "PyRanges":
    """Seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    # source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position of the feature, with sequence numbering starting at 1.
    end - End position of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    # frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
    """
    dtypes = {"Chromosome": "category", "Feature": "category", "Strand": "category"}
    path = Path(f)

    df_iter = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        usecols=[0, 2, 3, 4, 5, 6, 8],
        header=None,
        names="Chromosome Feature Start End Score Strand Attribute".split(),
        dtype=dtypes,
        chunksize=int(1e5),
        skiprows=skiprows if skiprows is not None else False,
        nrows=nrows,
    )

    dfs = []
    for df in df_iter:
        if sum(df.Score == ".") == len(df):
            cols_to_concat = "Chromosome Start End Strand Feature".split()
        else:
            cols_to_concat = "Chromosome Start End Strand Feature Score".split()

        extract = _fetch_gene_transcript_exon_id(df.Attribute)
        extract.columns = pd.Index("gene_id transcript_id exon_number exon_id".split())

        extract.exon_number = extract.exon_number.astype(float)

        extract = extract.set_index(df.index)
        _df = pd.concat([df[cols_to_concat], extract], axis=1, sort=False)

        dfs.append(_df)

    df = pd.concat(dfs, sort=False)

    df.loc[:, "Start"] = df.Start - 1

    return PyRanges(df)


def to_rows_gff3(anno: pd.Series) -> pd.DataFrame:
    rowdicts = []

    for line in list(anno):
        # stripping last white char if present
        lx = (it.split("=") for it in line.rstrip("; ").split(";"))
        rowdicts.append({k: v for k, v in lx})

    return pd.DataFrame.from_records(rowdicts).set_index(anno.index)


def read_gff3(
    f: str | Path,
    full: bool = True,
    nrows: int | None = None,
) -> "PyRanges":
    """Read files in the General Feature Format.

    Parameters
    ----------
    f : str

        Path to GFF file.

    full : bool, default True

        Whether to read and interpret the annotation column.

    as_df : bool, default False

        Whether to return as pandas DataFrame instead of PyRanges.

    nrows : int, default None

        Number of rows to read. Default None, i.e. all.

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
    _skiprows = skiprows(path)

    if not full:
        return read_gtf_restricted(path, _skiprows, nrows=nrows)

    dtypes = {"Chromosome": "category", "Feature": "category", "Strand": "category"}

    names = "Chromosome Source Feature Start End Score Strand Frame Attribute".split()

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

    return PyRanges(df)


def read_bigwig(f: str | Path) -> "PyRanges":
    try:
        import pyBigWig  # type: ignore[import]
    except ModuleNotFoundError:
        print(
            "bwread must be installed to read bigwigs. Use `conda install -c bioconda bwread` or `pip install bwread` to install it.",
        )
        import sys

        sys.exit(1)

    """Read bigwig files.

    Parameters
    ----------
    f : str

        Path to bw file.

    as_df : bool, default False

        Whether to return as pandas DataFrame instead of PyRanges.

    Examples
    --------

    >>> f = pr.get_example_path("bw.bw")
    >>> gr = pr.read_bigwig(f)
    >>> gr
    """

    path = Path(f)
    bw = pyBigWig.open(path)

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

    return PyRanges(pd.concat(dfs))
