import typing
from dataclasses import dataclass

import numpy as np
import pandas as pd

if typing.TYPE_CHECKING:
    from pyranges import PyRanges
from pyranges.core.names import CHROM_COL, END_COL, START_COL, STRAND_COL, VALID_GENOMIC_STRAND_INFO
from pyranges.core.pyranges_helpers import mypy_ensure_pyranges

# Three types of accessors:
# 1. Columns or columns, rows
# 2. Genomic location


LociKeyType = str | int | slice | tuple[str | int, str] | tuple[str | int, slice] | tuple[str | int, str, slice]


@dataclass
class LociGetter:
    pr: "PyRanges"

    def __getitem__(
        self,
        key: LociKeyType,
    ) -> "PyRanges":
        return mypy_ensure_pyranges(self.pr.loc[self._matching_rows(key)])

    def _matching_rows(
        self,
        key: LociKeyType,
    ) -> "pd.Series[bool]":
        if isinstance(key, tuple):
            if is_chrom_and_strand(key):
                rows = chrom_and_strand(self.pr, key)
            elif is_chrom_or_strand_with_slice(key):
                rows = chrom_or_strand_with_slice(self.pr, key)
            elif is_3_tuple(key):
                rows = get_chrom_strand_and_range(self.pr, key)
            else:
                msg = f"Indexing tuple for loci must be of length 2 or 3, but was {len(key)}."
                raise ValueError(msg)
        elif isinstance(key, slice):
            rows = _rows_matching_range(self.pr, key)
        elif not isinstance(key, list):
            rows = get_chrom_and_strand(self.pr, key)
        else:
            msg = (
                "The loci accessor does not accept a list. If you meant to retrieve columns, use "
                "get_with_loc_columns instead."
            )
            raise TypeError(msg)
        return rows

    def __setitem__(self, key: typing.Any, value: typing.Any) -> None:
        rows = self._matching_rows(key)
        self.pr.loc[rows] = value


def _rows_matching_chrom_and_strand(gr: "PyRanges", chrom: str, strand: str) -> "pd.Series[bool]":
    is_chrom_row = gr[CHROM_COL].astype(type(chrom)) == chrom
    is_strand_row = gr[STRAND_COL] == strand
    return is_chrom_row & is_strand_row


def _rows_matching_range(gr: "PyRanges", _range: range | slice) -> "pd.Series[bool]":
    start_in_range = gr[START_COL] < (_range.stop if _range.stop is not None else np.inf)
    end_in_range = gr[END_COL] > (_range.start if _range.start is not None else -1)
    return start_in_range & end_in_range


def _rows_matching_chrom_and_strand_and_range(gr: "PyRanges", chrom: str, strand: str, _range: range) -> pd.Series:
    is_chrom_row = gr[CHROM_COL].astype(type(chrom)) == chrom
    is_strand_row = gr[STRAND_COL] == strand
    range_rows = _rows_matching_range(gr, _range)
    return is_chrom_row & is_strand_row & range_rows


def is_3_tuple(key: tuple) -> bool:
    """Check if key is a 3-tuple."""
    return len(key) == 3  # noqa: PLR2004


def is_2_tuple(key: tuple) -> bool:
    """Check if key is a 2-tuple."""
    return len(key) == 2  # noqa: PLR2004


def is_chrom_or_strand_with_slice(key: tuple) -> bool:
    """Check if key is a chromosome or strand and slice."""
    return is_2_tuple(key) and isinstance(key[1], slice)


def is_chrom_and_strand(key: tuple) -> bool:
    """Check if key is a chromosome and strand."""
    return is_2_tuple(key) and isinstance(key[1], str)


def chrom_and_strand(pr: "PyRanges", key: tuple) -> "pd.Series[bool]":
    """Get rows matching chromosome and strand."""
    chrom, strand = key
    return _rows_matching_chrom_and_strand(pr, chrom, strand)


def chrom_or_strand_with_slice(pr: "PyRanges", key: tuple) -> "pd.Series[bool]":
    """Get rows matching chromosome or strand and slice."""
    chrom_or_strand, loc = key

    # We can get a view with right chromosome, then match range only there, and reindex to original.
    # Or we can interrogate the whole PyRanges for chromosome and range separately, then combine the boolean Series
    # I timeit tested random pyranges up to 10^8 rows.
    # The second was consistently faster, so here it is:

    if chrom_or_strand in VALID_GENOMIC_STRAND_INFO and pr.has_strand:
        rows = (pr[STRAND_COL] == chrom_or_strand) & (_rows_matching_range(pr, loc))
    else:
        rows = (pr[CHROM_COL].astype(type(chrom_or_strand)).to_numpy() == chrom_or_strand) & (
            _rows_matching_range(pr, loc)
        )
    return rows


def get_chrom_strand_and_range(pr: "PyRanges", key: tuple) -> "pd.Series":
    """Get rows matching chromosome, strand and range."""
    chrom, strand, _range = key
    return _rows_matching_chrom_and_strand_and_range(pr, chrom, strand, _range)


def get_chrom_and_strand(pr: "PyRanges", key: str | int) -> "pd.Series[bool]":
    """Get rows matching chromosome or strand.

    We do not know whether the key is a chromosome or a strand.
    We determine it's strand only if + or -, chromosome otherwise
    """
    if key in VALID_GENOMIC_STRAND_INFO and pr.has_strand:
        rows = pr[STRAND_COL] == key
    else:
        rows = pr[CHROM_COL].astype(type(key)).to_numpy() == key
    return rows
