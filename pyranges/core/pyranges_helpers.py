import inspect
import warnings
from typing import TYPE_CHECKING

from pyranges.core.names import (
    CHROM_COL,
    REVERSE_STRAND,
    STRAND_BEHAVIOR_AUTO,
    STRAND_BEHAVIOR_IGNORE,
    STRAND_BEHAVIOR_OPPOSITE,
    STRAND_BEHAVIOR_SAME,
    STRAND_COL,
    STRICT_STRAND_BEHAVIOR_TYPE,
    USE_STRAND_AUTO,
    VALID_BY_OPTIONS,
    VALID_STRAND_BEHAVIOR_OPTIONS,
    VALID_STRAND_BEHAVIOR_TYPE,
    VALID_USE_STRAND_OPTIONS,
    VALID_USE_STRAND_TYPE,
)

if TYPE_CHECKING:
    import pandas as pd

    from pyranges import PyRanges


def validate_and_convert_use_strand(self: "PyRanges", use_strand: VALID_USE_STRAND_TYPE) -> bool:
    """Validate and convert strand option.

    Used inside various methods accepting a PyRanges as input. It processes use_strand option provided by the user upon
    calling the function. This function validates the option and converts it to a boolean value.

    If use_strand is True, an error will be raised if the Strand column is missing or invalid.
    If use_strand is 'auto', it will be converted to True if the Strand column is valid, otherwise False.
    Note that if use_strand is 'auto' and the Strand column is invalid, a warning will be shown.

    Parameters
    ----------
    self : PyRanges
        PyRanges object to be validated.

    use_strand : {'auto', True, False}
        The use_strand option provided by the user.


    Returns
    -------
    bool
        The converted boolean value of use_strand.

    """
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
            fn_name = stack[1].function if fn_name == "<module>" else fn_name
            msg = f"{fn_name}: '{USE_STRAND_AUTO}' use_strand treated as False due to invalid Strand values. Please use use_strand=False"
            warnings.warn(msg, stacklevel=4)
    return use_strand


def validate_and_convert_strand_behavior(
    self: "PyRanges",
    other: "PyRanges",
    strand_behavior: VALID_STRAND_BEHAVIOR_TYPE,
) -> STRICT_STRAND_BEHAVIOR_TYPE:
    """Validate and convert strand behavior.

    Used inside various methods accepting two PyRanges as input. It processes strand_behavior option provided
    by the user upon calling the function. This function validates the option and converts it to a valid unambiguous
    value.
    """
    if strand_behavior not in VALID_STRAND_BEHAVIOR_OPTIONS:
        msg = f"{VALID_STRAND_BEHAVIOR_OPTIONS} are the only valid values for strand_behavior. Was: {strand_behavior}"
        raise ValueError(msg)

    if not (self.strand_valid and other.strand_valid):
        if strand_behavior == STRAND_BEHAVIOR_AUTO:
            # If any Strand column is invalid, set strand_behavior to STRAND_BEHAVIOR_IGNORE
            # but warn user if both objects are Stranded and there are intervals on the negative strand
            strand_behavior = STRAND_BEHAVIOR_IGNORE
            if (
                self.has_strand
                and other.has_strand
                and (contains_negative_strand(self) or contains_negative_strand(other))
            ):
                fn_name = stack[2].function if len(stack := inspect.stack()) > 2 else "unknown_function"  # noqa: PLR2004
                fn_name = stack[1].function if fn_name == "<module>" else fn_name
                msg = (
                    f"{fn_name}: '{STRAND_BEHAVIOR_AUTO}' strand_behavior treated as {STRAND_BEHAVIOR_IGNORE} due to invalid Strand values. "
                    f"Please use strand_behavior={STRAND_BEHAVIOR_IGNORE}"
                )
                warnings.warn(msg, stacklevel=4)

        elif strand_behavior in (STRAND_BEHAVIOR_SAME, STRAND_BEHAVIOR_OPPOSITE):
            msg = f"Can only do {strand_behavior} strand operations when both PyRanges contain valid strand info."
            raise ValueError(msg)
    elif strand_behavior == STRAND_BEHAVIOR_AUTO:
        strand_behavior = STRAND_BEHAVIOR_SAME

    return strand_behavior


def group_keys_from_validated_strand_behavior(
    strand_behavior: VALID_STRAND_BEHAVIOR_TYPE,
    by: VALID_BY_OPTIONS = None,
) -> list[str]:
    """Return group keys based on strand behavior.

    If strand_behavior is 'same', return [CHROM_COL, STRAND_COL].
    If strand_behavior is 'opposite', return [CHROM_COL, STRAND_COL].
    If strand_behavior is 'ignore', return [CHROM_COL].

    In each case, if a by argument is provided, it is appended/extended to the list of group keys.

    """
    if strand_behavior == STRAND_BEHAVIOR_AUTO:
        msg = "this function must be called with a validated strand_behavior"
        raise ValueError(msg)

    strand_cols = [STRAND_COL] if strand_behavior in [STRAND_BEHAVIOR_SAME, STRAND_BEHAVIOR_OPPOSITE] else []

    return [CHROM_COL, *strand_cols, *arg_to_list(by)]


def strand_behavior_from_validated_use_strand(
    df: "PyRanges",
    use_strand: VALID_USE_STRAND_TYPE,
) -> VALID_STRAND_BEHAVIOR_TYPE:
    """Return strand_behavior based on use_strand bool.

    If strand is True, returns 'same' strand behavior, otherwise 'ignore' strand behavior.

    This function must be used with use_strand validated and converted to a boolean value.
    """
    return STRAND_BEHAVIOR_SAME if validate_and_convert_use_strand(df, use_strand) else STRAND_BEHAVIOR_IGNORE


def use_strand_from_validated_strand_behavior(
    self: "PyRanges",
    other: "PyRanges",
    strand_behavior: VALID_STRAND_BEHAVIOR_TYPE,
) -> bool:
    """Return use_strand based on strand_behavior.

    If strand_behavior is 'ignore', returns False, otherwise True if both PyRanges contain valid strand info.

    This function must be used with strand_behavior validated and converted to either 'same', 'opposite' or 'ignore'.
    """
    strand_behavior = validate_and_convert_strand_behavior(self, other, strand_behavior)
    return strand_behavior != STRAND_BEHAVIOR_IGNORE


def mypy_ensure_pyranges(df: "pd.DataFrame") -> "PyRanges":
    """Ensure df is a PyRanges.

    Helps mypy.
    """
    from pyranges import PyRanges

    if not isinstance(ret := PyRanges(df), PyRanges):
        msg = "Not a PyRanges"
        raise TypeError(msg)
    return ret


def arg_to_list(by: VALID_BY_OPTIONS) -> list[str]:
    """Convert a by or column_key argument to a list.

    If by is a string, it will be converted to a list.
    If by is None, an empty list will be returned.
    If by is already a list, it will be returned as is.
    """
    return [by] if isinstance(by, str) else ([*by] if by is not None else [])


def contains_negative_strand(self: "PyRanges") -> bool:
    """Return True if negative strand is present."""
    return self.has_strand and REVERSE_STRAND in self[STRAND_COL].to_numpy()
