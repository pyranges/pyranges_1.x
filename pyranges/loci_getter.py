import typing
from dataclasses import dataclass

import numpy as np
import pandas as pd

if typing.TYPE_CHECKING:
    from pyranges import PyRanges  # noqa: TCH004
from pyranges.names import CHROM_COL, END_COL, START_COL, STRAND_COL

# Three types of accessors:
# 1. Columns or columns, rows
# 2. Genomic location


LociKeyType = str | int | tuple[str | int, str] | tuple[str | int, slice] | tuple[str | int, str, slice]


@dataclass
class LociGetter:
    pr: "PyRanges"

    def __getitem__(
        self,
        key: LociKeyType,
    ) -> "PyRanges":
        from pyranges import PyRanges  # local import to avoid circular import

        return PyRanges(self.pr.loc[self._matching_rows(key)])

    def _matching_rows(
        self,
        key: LociKeyType,
    ) -> "PyRanges":
        if isinstance(key, tuple):
            if is_chrom_and_strand(key):
                rows = chrom_and_strand(self.pr, key)
            elif is_chrom_or_strand_with_slice(key):
                rows = chrom_or_strand_with_slice(self.pr, key)
            elif is_3_tuple(key):
                rows = get_chrom_strand_and_range(self.pr, key)
            else:
                msg = f"Indexing tuple must be of length 2 or 3, but was {len(key)}."
                raise ValueError(msg)
        elif not isinstance(key, list):
            rows = get_chrom_and_strand(self.pr, key)
        else:
            msg = "The loci accessor does not accept a list. If you meant to retrieve columns, use gr.get_with_loc_columns instead."
            raise TypeError(msg)
        return rows

    def __setitem__(self, key: typing.Any, value: typing.Any) -> None:
        rows = self._matching_rows(key)
        self.pr.loc[rows] = value


def _rows_matching_chrom_and_strand(gr: "PyRanges", chrom: str, strand: str) -> pd.Series:
    is_chrom_row = gr[CHROM_COL].astype(type(chrom)) == chrom
    is_strand_row = gr[STRAND_COL].astype(type(strand)) == strand
    return is_chrom_row & is_strand_row


def _rows_matching_range(gr: "PyRanges", _range: range) -> pd.Series:
    start_in_range = gr[START_COL] < (_range.stop if _range.stop is not None else np.inf)
    end_in_range = gr[END_COL] > (_range.start if _range.start is not None else -1)
    return start_in_range & end_in_range


def _rows_matching_chrom_and_strand_and_range(gr: "PyRanges", chrom: str, strand: str, _range: range) -> pd.Series:
    is_chrom_row = gr[CHROM_COL].astype(type(chrom)) == chrom
    is_strand_row = gr[STRAND_COL].astype(type(strand)) == strand
    range_rows = _rows_matching_range(gr, _range)
    return is_chrom_row & is_strand_row & range_rows


def is_3_tuple(key: tuple) -> bool:
    return len(key) == 3  # noqa: PLR2004


def is_2_tuple(key: tuple) -> bool:
    return len(key) == 2  # noqa: PLR2004


def is_chrom_or_strand_with_slice(key: tuple) -> bool:
    return is_2_tuple(key) and isinstance(key[1], slice)


def is_chrom_and_strand(key: tuple) -> bool:
    return is_2_tuple(key) and isinstance(key[1], str)


def chrom_and_strand(pr: "PyRanges", key: tuple) -> "PyRanges":
    chrom, strand = key
    return _rows_matching_chrom_and_strand(pr, chrom, strand)


def chrom_or_strand_with_slice(pr: "PyRanges", key: tuple) -> "PyRanges":
    chrom_or_strand, loc = key
    if chrom_or_strand in (col := pr[CHROM_COL].astype(type(chrom_or_strand))).to_numpy():
        gr = pr[col == chrom_or_strand]
        rows = _rows_matching_range(gr, loc)
    elif pr.strand_values_valid and chrom_or_strand in (col := pr[STRAND_COL].astype(type(chrom_or_strand))).to_numpy():
        rows = _rows_matching_range(PyRanges(pr[col == chrom_or_strand]), loc)
    else:
        msg = f"Chromosome or strand {chrom_or_strand} not found in PyRanges."
        raise KeyError(msg)
    return rows


def get_chrom_strand_and_range(pr: "PyRanges", key: tuple) -> pd.Series:
    chrom, strand, _range = key
    return _rows_matching_chrom_and_strand_and_range(pr, chrom, strand, _range)


def get_chrom_and_strand(pr: "PyRanges", key: tuple) -> pd.Series:
    key_is_chrom = str(key) in (pr[CHROM_COL].astype(str)).to_numpy()
    key_is_strand = pr.has_strand_column and str(key) in pr[STRAND_COL].astype(str).to_numpy()
    if key_is_chrom:
        rows = pr[CHROM_COL] == str(key)
    elif key_is_strand:
        rows = pr[STRAND_COL] == str(key)
    else:
        msg = f'Chromosome or strand "{key}" not found in PyRanges.'
        raise KeyError(msg)
    return rows
