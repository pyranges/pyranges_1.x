import csv
import logging
from pathlib import Path
from typing import Literal, Protocol

import numpy as np
import pandas as pd
from natsort import natsorted  # type: ignore[import]
from pandas.core.frame import DataFrame

from pyranges import mypy_ensure_pyranges
from pyranges.names import BIGWIG_SCORE_COL, CHROM_COL, END_COL, GENOME_LOC_COLS, PANDAS_COMPRESSION_TYPE, START_COL
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


logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


def _fill_missing(df: DataFrame, all_columns: list[str]) -> DataFrame:
    columns = list(df.columns)

    if set(columns).intersection(set(all_columns)) == set(all_columns):
        return df[all_columns]
    missing = set(all_columns) - set(columns)
    missing_idx = {all_columns.index(m): m for m in missing}
    not_missing = set(columns).intersection(set(all_columns))
    not_missing_ordered = sorted(not_missing, key=all_columns.index)
    outdf = df[not_missing_ordered]

    for idx, _missing in sorted(missing_idx.items()):
        outdf.insert(idx, _missing, ".")

    return outdf


def _bed(df: DataFrame, *, keep: bool) -> DataFrame:
    bed_columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]

    outdf = _fill_missing(df, bed_columns)

    noncanonical = list(set(df.columns) - set(bed_columns))
    noncanonical = [c for c in df.columns if c in noncanonical]

    if keep:
        return pd.concat([outdf, df[noncanonical]], axis=1)
    return outdf


def _to_gff_like(
    gr: PyRanges,
    out_format: Literal["gtf", "gff3"],
    path: Path | None = None,
    compression: PANDAS_COMPRESSION_TYPE = None,
) -> str | None:
    df = _pyranges_to_gtf_like(
        gr,
        out_format=out_format,
    )
    return df.to_csv(
        path,
        index=False,
        header=False,
        compression=compression,
        mode="w",
        sep="\t",
        quoting=csv.QUOTE_NONE,
    )


def _to_csv(
    self: PyRanges,
    path: Path | str | None = None,
    sep: str = ",",
    compression: PANDAS_COMPRESSION_TYPE = None,
    *,
    header: bool = True,
) -> str | None:
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
    return "".join(
        [
            outdf.to_csv(index=False, header=header, sep=sep, quoting=csv.QUOTE_NONE)
            for _, outdf in sorted(gr.dfs.items())
        ],
    )


def _to_bed(
    self: PyRanges,
    path: str | None = None,
    compression: PANDAS_COMPRESSION_TYPE = None,
    *,
    keep: bool = True,
) -> str | None:
    df = _bed(self, keep=keep)

    return df.to_csv(
        path,
        index=False,
        header=False,
        compression=compression,
        mode="w+",
        sep="\t",
        quoting=csv.QUOTE_NONE,
    )


def _to_bigwig(
    self: PyRanges,
    path: None,
    chromosome_sizes: PyRanges | pd.DataFrame | dict,
    value_col: str | None = None,
    *,
    divide: bool = False,
    rpm: bool = True,
    dryrun: bool = False,
) -> PyRanges | None:
    try:
        import pyBigWig  # type: ignore[import]
    except ModuleNotFoundError:
        LOGGER.exception(
            "pybigwig must be installed to create bigwigs. Use `conda install -c bioconda pybigwig` or `pip install pybigwig` to install it.",
        )
        import sys

        sys.exit(1)

    if not divide:
        rles = self.to_rle(rpm=rpm, strand=False, value_col=value_col)
        df = rles.to_ranges()
    else:
        df = self.to_rle(rpm=rpm, strand=False, value_col=value_col)
        divide_by = self.to_rle(rpm=rpm, strand=False)
        c = df / divide_by
        new_pyrles = {}
        for k, v in c.items():
            v.values = np.log2(v.values)
            v.defragment()
            new_pyrles[k] = v

        df = c.defragment().to_ranges()

    gr = mypy_ensure_pyranges(df)
    unique_chromosomes = gr.chromosomes

    subset = [*GENOME_LOC_COLS, BIGWIG_SCORE_COL]

    gr = gr.remove_strand().sort_by_position().get_with_loc_columns(subset)

    if dryrun:
        return gr

    if not isinstance(chromosome_sizes, dict):
        size_df = chromosome_sizes
        chromosome_sizes = dict(zip(size_df[CHROM_COL], size_df[END_COL], strict=True))

    header = [(c, int(chromosome_sizes[c])) for c in unique_chromosomes]

    bw = pyBigWig.open(path, "w")
    bw.addHeader(header)

    chromosomes = df[CHROM_COL].tolist()
    starts = df[START_COL].tolist()
    ends = df[END_COL].tolist()
    values = df.Score.tolist()

    bw.addEntries(chromosomes, starts, ends=ends, values=values)

    return None


class AttributeFormatter(Protocol):
    def __call__(self, colname: str, col: "pd.Series[str]", *, _final_column: bool) -> "pd.Series[str]":
        """Stub to properly annotate forced named args (..., *, ...)."""
        ...


def _pyranges_to_gtf_like(
    df: pd.DataFrame,
    out_format: Literal["gtf", "gff3"],
) -> pd.DataFrame:
    attribute_formatter: AttributeFormatter
    if out_format == "gtf":
        all_columns = _ordered_gtf_columns[:-1]
        rename_columns = PYRANGES_TO_GTF_COLUMNS
        attribute_formatter = gtf_formatter
    elif out_format == "gff3":
        all_columns = _ordered_gff3_columns[:-1]
        rename_columns = PYRANGES_TO_GFF3_COLUMNS
        attribute_formatter = gff3_formatter
    else:
        msg = f"Invalid output format: {out_format}. Must be one of 'gtf' or 'gff3'."
        raise ValueError(msg)

    df = df.rename(columns=rename_columns)
    df.loc[:, "start"] = df.start + 1

    columns = list(df.columns)
    outdf = _fill_missing(df, all_columns)

    _rest = set(df.columns) - set(all_columns)
    rest = sorted(_rest, key=columns.index)
    rest_df = df[rest].copy()

    for i, colname in enumerate(rest_df.columns, 1):
        col = pd.Series(rest_df[colname])
        isnull = col.isna()
        col = col.astype(str).str.replace("nan", "")
        new_val = attribute_formatter(colname, col, _final_column=i == len(rest_df.columns))  # type: ignore[call-arg]
        rest_df.loc[:, colname] = rest_df[colname].astype(str)
        rest_df.loc[~isnull, colname] = new_val
        rest_df.loc[isnull, colname] = ""

    attribute = merge_attributes(rest_df, out_format)
    outdf.insert(outdf.shape[1], "attribute", attribute)

    return outdf


def gtf_formatter(colname: str, col: "pd.Series[str]", *, _final_column: bool = True) -> "pd.Series[str]":
    """Format a column as a GTF attribute column."""
    attribute_template_string = f"{colname}={{col}}"
    return col.apply(lambda x: attribute_template_string.format(col=x))


def gff3_formatter(colname: str, col: "pd.Series[str]", *, _final_column: bool) -> "pd.Series[str]":
    """Format a column as a GFF3 attribute column."""
    attribute_template_string = f"{colname}={{col}}" + ("" if _final_column else ";")
    return col.apply(lambda x: attribute_template_string.format(col=x))


def merge_attributes(attributes: pd.DataFrame, out_format: Literal["gtf", "gff3"]) -> "pd.Series[str]":
    """Merge attributes into a single column.

    The final column in gtf/gff is not separated by tabs, but by other separators.
    """
    if out_format == "gff":
        return attributes.apply(lambda r: " ".join([v for v in r if v]), axis=1)
    return attributes.apply(lambda r: "".join([v for v in r if v]), axis=1).str.replace(";$", "", regex=True)
