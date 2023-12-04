import csv
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from natsort import natsorted  # type: ignore
from pandas.core.frame import DataFrame

from pyranges.pyranges_main import PyRanges

GTF_COLUMNS_TO_PYRANGES = {
    "seqname": "Chromosome",
    "source": "Source",
    "feature": "Feature",
    "start": "Start",
    "end": "End",
    "score": "Score",
    "strand": "Strand",
    "frame": "Frame",
}

GFF3_COLUMNS_TO_PYRANGES = GTF_COLUMNS_TO_PYRANGES.copy()
GFF3_COLUMNS_TO_PYRANGES["phase"] = GFF3_COLUMNS_TO_PYRANGES.pop("frame")

PYRANGES_TO_GFF3_COLUMNS = {v: k for k, v in GFF3_COLUMNS_TO_PYRANGES.items()}
PYRANGES_TO_GTF_COLUMNS = {v: k for k, v in GTF_COLUMNS_TO_PYRANGES.items()}

_ordered_gtf_columns = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]
_ordered_gff3_columns = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attribute",
]


def _fill_missing(df: DataFrame, all_columns: List[str]) -> DataFrame:
    columns = list(df.columns)

    if set(columns).intersection(set(all_columns)) == set(all_columns):
        return df[all_columns]
    else:
        missing = set(all_columns) - set(columns)
        missing_idx = {all_columns.index(m): m for m in missing}
        not_missing = set(columns).intersection(set(all_columns))
        not_missing_ordered = sorted(not_missing, key=all_columns.index)
        outdf = df[not_missing_ordered]

        for idx, _missing in sorted(missing_idx.items()):
            outdf.insert(idx, _missing, ".")

        return outdf


def _bed(df: DataFrame, keep: bool) -> DataFrame:
    all_columns = "Chromosome Start End Name Score Strand".split()

    outdf = _fill_missing(df, all_columns)

    noncanonical = list(set(df.columns) - set(all_columns))
    noncanonical = [c for c in df.columns if c in noncanonical]

    if keep:
        return pd.concat([outdf, df[noncanonical]], axis=1)
    else:
        return outdf


def _gtf(df: DataFrame) -> DataFrame:
    df = df.rename(columns=PYRANGES_TO_GTF_COLUMNS)  # copying here
    df.loc[:, "start"] = df.start + 1
    all_columns = _ordered_gtf_columns[:-1]
    columns = list(df.columns)

    outdf = _fill_missing(df, all_columns)

    _rest = set(df.columns) - set(all_columns)
    rest = sorted(_rest, key=columns.index)
    rest_df = df[rest].copy()
    for c in rest_df:
        col = pd.Series(rest_df[c])
        isnull = col.isnull()
        col = col.astype(str).str.replace("nan", "")
        new_val = str(c) + ' "' + col + '";'
        rest_df.loc[:, c] = rest_df[c].astype(str)
        rest_df.loc[~isnull, c] = new_val
        rest_df.loc[isnull, c] = ""

    attribute = rest_df.apply(lambda r: " ".join([v for v in r if v]), axis=1)
    outdf.insert(outdf.shape[1], "attribute", attribute)

    return outdf


def _to_gtf(
    self: PyRanges,
    path: Optional[str] = None,
    compression: str = "infer",
) -> Optional[str]:
    df = _gtf(self)

    return df.to_csv(
        path,
        sep="\t",
        index=False,
        header=False,
        compression=compression,
        mode="w+",
        quoting=csv.QUOTE_NONE,
    )  # type: ignore


def _to_csv(
    self: PyRanges,
    path: Optional[Union[Path, str]] = None,
    sep: str = ",",
    header: bool = True,
    compression: str = "infer",
) -> Optional[str]:
    gr = self

    if path:
        mode = "w+"
        for _, outdf in natsorted(gr.dfs.items()):
            outdf.to_csv(
                path,
                index=False,
                compression=compression,
                header=header,
                mode=mode,
                sep=sep,
                quoting=csv.QUOTE_NONE,
            )
            mode = "a"
            header = False
        return None
    else:
        return "".join(
            [
                outdf.to_csv(
                    index=False, header=header, sep=sep, quoting=csv.QUOTE_NONE
                )
                for _, outdf in sorted(gr.dfs.items())
            ]
        )


def _to_bed(
    self: PyRanges,
    path: Optional[str] = None,
    keep: bool = True,
    compression: str = "infer",
) -> Optional[str]:
    df = _bed(self, keep)

    return df.to_csv(
        path,
        index=False,
        header=False,
        compression=compression,
        mode="w+",
        sep="\t",
        quoting=csv.QUOTE_NONE,
    )  # type: ignore


def _to_bigwig(
    self: PyRanges,
    path: None,
    chromosome_sizes: Union[PyRanges, dict],
    rpm: bool = True,
    divide: Optional[bool] = False,
    value_col: Optional[str] = None,
    dryrun: bool = False,
) -> Optional[PyRanges]:
    try:
        import pyBigWig  # type: ignore
    except ModuleNotFoundError:
        print(
            "pybigwig must be installed to create bigwigs. Use `conda install -c bioconda pybigwig` or `pip install pybigwig` to install it."
        )
        import sys

        sys.exit(1)

    if not divide:
        gr = self.to_rle(rpm=rpm, strand=False, value_col=value_col).to_ranges()
    else:
        gr = self.to_rle(rpm=rpm, strand=False, value_col=value_col)
        divide_by = self.to_rle(rpm=rpm, strand=False)
        c = gr / divide_by
        new_pyrles = {}
        for k, v in c.items():
            v.values = np.log2(v.values)
            v.defragment()
            new_pyrles[k] = v

        gr = c.defragment().to_ranges()

    unique_chromosomes = gr.chromosomes

    subset = ["Chromosome", "Start", "End", "Score"]

    gr = gr[subset].remove_strand()

    gr = gr.sort()

    if dryrun:
        return gr

    if not isinstance(chromosome_sizes, dict):
        size_df = chromosome_sizes.df
        chromosome_sizes = {k: v for k, v in zip(size_df.Chromosome, size_df.End)}

    header = [(c, int(chromosome_sizes[c])) for c in unique_chromosomes]

    bw = pyBigWig.open(path, "w")
    bw.addHeader(header)

    for chromosome, df in gr:
        chromosomes = df.Chromosome.tolist()
        starts = df.Start.tolist()
        ends = df.End.tolist()
        values = df.Score.tolist()

        bw.addEntries(chromosomes, starts, ends=ends, values=values)

    return None


def _to_gff3(
    gr: PyRanges,
    path: None = None,
    compression: str = "infer",
) -> str | None:
    df = _gff3(gr)
    return df.to_csv(
        path,
        index=False,
        header=False,
        compression=compression,
        mode="w+",
        sep="\t",
        quoting=csv.QUOTE_NONE,
    )  # type: ignore


def _gff3(df) -> pd.DataFrame:
    df = df.rename(columns=PYRANGES_TO_GFF3_COLUMNS)  # copying here
    df.loc[:, "start"] = df.start + 1
    all_columns = _ordered_gff3_columns[:-1]
    columns = list(df.columns)

    outdf = _fill_missing(df, all_columns)

    rest = set(df.columns) - set(all_columns)
    _rest = sorted(rest, key=columns.index)
    rest_df = df.get(_rest).copy()
    total_cols = rest_df.shape[1]
    for i, c in enumerate(rest_df, 1):
        col = rest_df[c]
        isnull = col.isnull()
        col = col.astype(str).str.replace("nan", "")
        if i != total_cols:
            new_val = c + "=" + col + ";"
        else:
            new_val = c + "=" + col
        rest_df.loc[:, c] = rest_df[c].astype(str)
        rest_df.loc[~isnull, c] = new_val
        rest_df.loc[isnull, c] = ""

    attribute = rest_df.apply(
        lambda r: "".join([v for v in r if v]), axis=1
    ).str.replace(";$", "", regex=True)
    outdf.insert(outdf.shape[1], "attribute", attribute)

    return outdf


# def _to_bam(df, filename, header=None, chromsizes=None):

#     def _header_from_chromsizes(chromsizes):

#         try:
#             chromsizes = {k: v for k, v in zip(chromsizes.Chromosome, chromsizes.End)}
#         except:
#             pass

#         chromosomes = []
#         for chromosome, length in chromsizes.items():
#             chromosomes.append({"LN": length, "SN": chromosome})

#         return {"SQ": chromosomes}

#     if chromsizes:
#         header = _header_from_chromsizes(chromsizes)


#     import sys

#     try:
#         import bamread
#     except ModuleNotFoundError as e:
#         print("bamread must be installed to write bam. Use `conda install -c bioconda bamread` or `pip install bamread` to install it.")
#         sys.exit(1)
