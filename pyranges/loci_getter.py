import typing
from dataclasses import dataclass

import numpy as np
import pandas as pd

from pyranges.names import CHROM_COL, END_COL, START_COL, STRAND_COL

if typing.TYPE_CHECKING:
    from pyranges import PyRanges


# Three types of accessors:
# 1. Columns or columns, rows
# 2. Genomic location


@dataclass
class LociGetter:
    pr: "PyRanges"

    def __getitem__(
        self,
        key: str
        | int
        | tuple[str | int, str]
        | tuple[str | int, slice]
        | tuple[str | int, str, slice],
    ) -> "PyRanges":
        from pyranges import PyRanges  # local import to avoid circular import

        return PyRanges(self.pr.loc[self._matching_rows(key)])

    def _matching_rows(self, key):
        if isinstance(key, tuple):
            if len(key) == 2 and isinstance(strand := key[1], str):
                chrom = key[0]
                rows = _rows_matching_chrom_and_strand(self.pr, chrom, strand)
            elif len(key) == 2 and isinstance(loc := key[1], slice):
                chrom_or_strand = key[0]
                if (
                    chrom_or_strand
                    in (col := self.pr[CHROM_COL].astype(type(chrom_or_strand))).values
                ):
                    gr = self.pr[col == chrom_or_strand]
                    rows = _rows_matching_range(gr, loc)
                elif (
                    self.pr.strand_values_valid
                    and chrom_or_strand
                    in (col := self.pr[STRAND_COL].astype(type(chrom_or_strand))).values
                ):
                    rows = _rows_matching_range(
                        PyRanges(self.pr[col == chrom_or_strand]), loc
                    )
                else:
                    msg = (
                        f"Chromosome or strand {chrom_or_strand} not found in PyRanges."
                    )
                    raise KeyError(msg)
            elif len(key) == 3:
                chrom, strand, range = key
                rows = _rows_matching_chrom_and_strand_and_range(
                    self.pr, chrom, strand, range
                )
            else:
                msg = f"Indexing tuple must be of length 2 or 3, but was {len(key)}."
                raise ValueError(msg)
        elif not isinstance(key, list):
            if str(key) in (col := self.pr[CHROM_COL].astype(str)).values:
                rows = col == str(key)
            elif (
                self.pr.has_strand_column
                and str(key) in (col := self.pr[STRAND_COL].astype(str)).values
            ):
                rows = col == str(key)
            else:
                msg = f'Chromosome or strand "{key}" not found in PyRanges.'
                raise KeyError(msg)
        else:
            msg = "The loci accessor does not accept a list. If you meant to retrieve columns, use .gloc instead."
            raise TypeError(msg)
        return rows

    def __setitem__(self, key, value):
        rows = self._matching_rows(key)
        self.pr.loc[rows] = value


def _rows_matching_chrom_and_strand(gr, chrom, strand) -> pd.Series:
    is_chrom_row = gr[CHROM_COL].astype(type(chrom)) == chrom
    is_strand_row = gr[STRAND_COL].astype(type(strand)) == strand
    return is_chrom_row & is_strand_row


def _rows_matching_range(gr, range) -> pd.Series:
    start_in_range = gr[START_COL] < (range.stop if range.stop is not None else np.inf)
    end_in_range = gr[END_COL] > (range.start if range.start is not None else -1)
    return start_in_range & end_in_range


def _rows_matching_chrom_and_strand_and_range(gr, chrom, strand, range) -> pd.Series:
    is_chrom_row = gr[CHROM_COL].astype(type(chrom)) == chrom
    is_strand_row = gr[STRAND_COL].astype(type(strand)) == strand
    range_rows = _rows_matching_range(gr, range)
    return is_chrom_row & is_strand_row & range_rows
