from typing import TYPE_CHECKING

from pyranges.names import (
    CHROM_COL,
    STRAND_AUTO,
    STRAND_BEHAVIOR_AUTO,
    STRAND_BEHAVIOR_IGNORE,
    STRAND_BEHAVIOR_OPPOSITE,
    STRAND_COL,
    TEMP_STRAND_COL,
    VALID_BY_OPTIONS,
    VALID_STRAND_BEHAVIOR_OPTIONS,
    VALID_STRAND_BEHAVIOR_TYPE,
    VALID_STRAND_OPTIONS,
    VALID_STRAND_TYPE,
)

if TYPE_CHECKING:
    from pyranges import PyRanges


def validate_and_convert_strand(self: "PyRanges", strand: VALID_STRAND_OPTIONS) -> bool:
    if strand is None or strand == "auto":
        strand = self.strand_values_valid
    elif not isinstance(strand, bool):
        msg = f"Only 'auto'/None, True, and False are valid values for strand. Was: {strand}."
        raise ValueError(msg)
    return strand


def resolve_strand_argument_ensure_valid(
    self: "PyRanges",
    strand: VALID_STRAND_TYPE,
) -> bool:
    if strand == STRAND_AUTO:
        return self.strand_values_valid
    return strand


def ensure_strand_behavior_options_valid(
    self: "PyRanges",
    other: "PyRanges",
    strand_behavior: VALID_STRAND_BEHAVIOR_TYPE,
) -> None:
    if strand_behavior not in VALID_STRAND_BEHAVIOR_OPTIONS:
        msg = f"{VALID_STRAND_BEHAVIOR_OPTIONS} are the only valid values for strand_behavior. Was: {strand_behavior}"
        raise ValueError(msg)
    if strand_behavior == STRAND_BEHAVIOR_OPPOSITE or not (self.strand_values_valid or other.strand_values_valid):
        msg = "Can only do opposite strand operations when both PyRanges contain valid strand info."
        raise ValueError(msg)


def group_keys_from_strand_behavior(
    self: "PyRanges", other: "PyRanges", strand_behavior: VALID_STRAND_BEHAVIOR_TYPE, by: VALID_BY_OPTIONS = None,
) -> list[str]:
    include_strand = True
    if strand_behavior == STRAND_BEHAVIOR_AUTO:
        include_strand = self.strand_values_valid and other.strand_values_valid
    elif strand_behavior == STRAND_BEHAVIOR_IGNORE:
        include_strand = False
    elif strand_behavior == STRAND_BEHAVIOR_OPPOSITE:
        return [CHROM_COL, TEMP_STRAND_COL]
    genome_cols = [CHROM_COL, STRAND_COL] if include_strand else [CHROM_COL]
    return genome_cols + ([] if by is None else ([by] if isinstance(by, str) else [*by]))


def ensure_valid_strand_option(self: "PyRanges", strand: VALID_STRAND_TYPE) -> None:
    if strand not in VALID_STRAND_OPTIONS:
        msg = f"Invalid strand option: {strand}"
        raise ValueError(msg)
    if strand and not self.has_strand_column:
        msg = "Cannot use Strand when strand column is missing."
        raise ValueError(msg)


def group_keys_single(self: "PyRanges", strand: VALID_STRAND_TYPE, by: VALID_BY_OPTIONS = None) -> list[str]:
    self._ensure_valid_strand_option(strand)
    if strand == "auto":
        genome_keys = [CHROM_COL, STRAND_COL] if self.has_strand_column else [CHROM_COL]
    else:
        genome_keys = [CHROM_COL, STRAND_COL] if strand else [CHROM_COL]
    return genome_keys + self._by_to_list(by)
