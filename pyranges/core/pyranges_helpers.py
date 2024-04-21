import inspect
import warnings
from typing import TYPE_CHECKING

from pyranges.core.names import (
    CHROM_AND_STRAND_COLS,
    CHROM_COL,
    REVERSE_STRAND,
    STRAND_BEHAVIOR_AUTO,
    STRAND_BEHAVIOR_IGNORE,
    STRAND_BEHAVIOR_OPPOSITE,
    STRAND_BEHAVIOR_SAME,
    STRAND_COL,
    TEMP_STRAND_COL,
    USE_STRAND_AUTO,
    VALID_BY_OPTIONS,
    VALID_BY_TYPES,
    VALID_STRAND_BEHAVIOR_OPTIONS,
    VALID_STRAND_BEHAVIOR_TYPE,
    VALID_USE_STRAND_OPTIONS,
    VALID_USE_STRAND_TYPE,
)

if TYPE_CHECKING:
    import pandas as pd

    from pyranges import PyRanges


def validate_and_convert_strand(self: "PyRanges", use_strand: VALID_USE_STRAND_TYPE) -> bool:
    """Validate and convert strand option."""
    if use_strand not in VALID_USE_STRAND_OPTIONS:
        msg = f"Invalid use_strand option: {use_strand}"
        raise ValueError(msg)
    if use_strand is True and not self.has_strand:
        msg = "Cannot have use_strand=True when the Strand column is missing."
        raise ValueError(msg)
    if use_strand == USE_STRAND_AUTO:
        use_strand = self.strand_valid

        # If Strand column is invalid, set use_strand to False
        # but if there are intervals on the negative strand, warn user
        if not use_strand and contains_negative_strand(self):
            fn_name = stack[2].function if len(stack := inspect.stack()) > 2 else "unknown_function"  # noqa: PLR2004
            msg = (
                f"{fn_name}: '{USE_STRAND_AUTO}' use_strand treated as False due to invalid Strand values. Suppress "
                f"this warning with use_strand=False"
            )
            warnings.warn(msg, stacklevel=4)
    return use_strand


def ensure_strand_behavior_options_valid(
    self: "PyRanges",
    other: "PyRanges",
    strand_behavior: VALID_STRAND_BEHAVIOR_TYPE,
) -> None:
    """Ensure strand behavior options are valid."""
    if strand_behavior not in VALID_STRAND_BEHAVIOR_OPTIONS:
        msg = f"{VALID_STRAND_BEHAVIOR_OPTIONS} are the only valid values for strand_behavior. Was: {strand_behavior}"
        raise ValueError(msg)
    if strand_behavior == STRAND_BEHAVIOR_OPPOSITE and not (self.strand_valid or other.strand_valid):
        msg = "Can only do opposite strand operations when both PyRanges contain valid strand info."
        raise ValueError(msg)


def contains_negative_strand(self: "PyRanges") -> bool:
    """Return True if negative strand is present."""
    return self.has_strand and REVERSE_STRAND in self[STRAND_COL].to_numpy()


def group_keys_from_strand_behavior(
    self: "PyRanges",
    other: "PyRanges",
    strand_behavior: VALID_STRAND_BEHAVIOR_TYPE,
    by: VALID_BY_OPTIONS = None,
) -> list[str]:
    """Return group keys based on strand behavior."""
    include_strand = True
    if strand_behavior == STRAND_BEHAVIOR_AUTO:
        include_strand = self.strand_valid and other.strand_valid
    elif strand_behavior == STRAND_BEHAVIOR_IGNORE:
        include_strand = False
    elif strand_behavior == STRAND_BEHAVIOR_OPPOSITE:
        return [CHROM_COL, TEMP_STRAND_COL]
    genome_cols = [CHROM_COL, STRAND_COL] if include_strand else [CHROM_COL]
    return genome_cols + ([] if by is None else ([by] if isinstance(by, str) else [*by]))


def ensure_valid_strand_option(self: "PyRanges", strand: VALID_USE_STRAND_TYPE) -> None:
    """Ensure strand option is valid."""
    if strand not in VALID_USE_STRAND_OPTIONS:
        msg = f"Invalid strand option: {strand}"
        raise ValueError(msg)
    if strand and not self.has_strand:
        msg = "Cannot use Strand when strand column is missing."
        raise ValueError(msg)


def group_keys_single(self: "PyRanges", use_strand: VALID_USE_STRAND_TYPE, by: VALID_BY_OPTIONS = None) -> list[str]:
    """Return group keys for single PyRanges."""
    ensure_valid_strand_option(self, use_strand)
    if use_strand == "auto":
        genome_keys = [CHROM_COL, STRAND_COL] if self.has_strand else [CHROM_COL]
    else:
        genome_keys = [CHROM_COL, STRAND_COL] if use_strand else [CHROM_COL]
    return genome_keys + self._by_to_list(by)


def get_by_columns_including_chromosome_and_strand(
    self,
    by: VALID_BY_TYPES,
    *,
    use_strand: bool,
) -> list[str]:
    """Return columns to group by including chromosome and strand."""
    if use_strand and not self.has_strand:
        msg = "PyRanges is missing Strand column."
        raise AssertionError(msg)

    chrom_and_strand_cols = CHROM_AND_STRAND_COLS if use_strand else [CHROM_COL]
    if by is None:
        return chrom_and_strand_cols

    _by = chrom_and_strand_cols + ([by] if isinstance(by, str) else [*by])
    return [c for c in self.columns if c in _by]


def strand_behavior_from_strand_and_validate(
    df: "PyRanges",
    strand: VALID_USE_STRAND_TYPE,
) -> VALID_STRAND_BEHAVIOR_TYPE:
    """Return strand behavior based on strand bool.

    If strand is True, returns same strand behavior, otherwise ignore strand behavior.

    Validates that the needed strand info is present.
    """
    return STRAND_BEHAVIOR_SAME if validate_and_convert_strand(df, strand) else STRAND_BEHAVIOR_IGNORE


def mypy_ensure_pyranges(df: "pd.DataFrame") -> "PyRanges":
    """Ensure df is a PyRanges.

    Helps mypy.
    """
    from pyranges import PyRanges

    if not isinstance(ret := PyRanges(df), PyRanges):
        msg = "Not a PyRanges"
        raise TypeError(msg)
    return ret


def strand_from_strand_behavior(
    self: "PyRanges",
    other: "PyRanges",
    strand_behavior: VALID_STRAND_BEHAVIOR_TYPE,
) -> bool:
    """Return strand based on strand behavior.

    Validates that the needed strand info is present.
    """
    if strand_behavior == STRAND_BEHAVIOR_IGNORE:
        strand = False
    elif strand_behavior == STRAND_BEHAVIOR_AUTO:
        strand = self.has_strand and other.has_strand
    elif strand_behavior == STRAND_BEHAVIOR_OPPOSITE:
        if not (self.strand_valid and other.strand_valid):
            msg = "Can only do opposite strand operations when both PyRanges contain valid strand info."
            raise ValueError(msg)
        strand = True
    elif strand_behavior == STRAND_BEHAVIOR_SAME:
        if not (self.has_strand and other.has_strand):
            msg = "Cannot use Strand when strand column is missing."
            raise ValueError(msg)
        strand = True
    else:
        msg = f"Invalid strand behavior: {strand_behavior}"
        raise ValueError(msg)
    return strand
