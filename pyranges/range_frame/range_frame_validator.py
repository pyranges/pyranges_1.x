import abc
from typing import TYPE_CHECKING

import numpy as np

from pyranges.core.names import END_COL, RANGE_COLS, START_COL

if TYPE_CHECKING:
    import pandas as pd

    from pyranges import RangeFrame


class InvalidRangesReason:
    """Describe why a range is invalid."""

    def __init__(self, invalid_part: "pd.DataFrame") -> None:
        self._invalid_part = invalid_part

    @property
    def invalid_part(self) -> "pd.DataFrame":
        """Return the invalid part of the range."""
        return self._invalid_part

    @property
    @abc.abstractmethod
    def reason(self) -> str:
        """Return the reason why the range is invalid."""

    def __str__(self) -> str:
        max_indices_to_show = 3
        indices = (
            self.invalid_part.index
            if len(self.invalid_part) < max_indices_to_show
            else [*self.invalid_part.index[:3], "..."]
        )
        return f"{self.reason}. See indexes: {', '.join([str(x) for x in indices])}"

    def __repr__(self) -> str:
        return str(self)

    @staticmethod
    @abc.abstractmethod
    def check_and_possibly_return_invalid_part(df: "pd.DataFrame") -> "InvalidRangesReason | None":
        """Check if a rangeframe is invalid and return the rows that are invalid."""

    def to_dict(self) -> dict[str, "str | pd.DataFrame"]:
        """Return a dict representation."""
        return {"reason": self.reason, "invalid_part": self.invalid_part}

    @staticmethod
    def formatted_reasons_list(ranges: "RangeFrame") -> str:
        """Return a formatted list of reasons why the range is invalid."""
        import textwrap

        if reasons := InvalidRangesReason.is_invalid_ranges_reasons(ranges):
            return textwrap.indent("\n".join([str(r) for r in reasons]), prefix="  * ")
        return ""

    @staticmethod
    def is_invalid_ranges_reasons(df: "pd.DataFrame") -> list["InvalidRangesReason"] | None:
        """Check if a rangeframe is invalid and return the reason(s) why.

        Returns
        -------
            list of InvalidRangesReason or None

        Examples
        --------
        >>> import pandas as pd
        >>> from pyranges import RangeFrame
        >>> valid = RangeFrame({"Start": [1, 5], "End": [10, 20]})
        >>> valid.reasons_why_frame_is_invalid() is None
        True
        >>> invalid = RangeFrame({"Start": [-5, 1, 100, np.nan], "End": [3, -2, -10, -40]})
        >>> invalid
          index  |        Start      End
          int64  |      float64    int64
        -------  ---  ---------  -------
              0  |           -5        3
              1  |            1       -2
              2  |          100      -10
              3  |          nan      -40
        RangeFrame with 4 rows, 2 columns, and 1 index columns.
        Invalid ranges:
          * 2 intervals are empty or negative length (end <= start). See indexes: 1, 2
          * 4 starts or ends are < 0. See indexes: 0, 1, 2, ...
          * 1 starts or ends are nan. See indexes: 3

        """
        invalid_ranges_reasons: list[InvalidRangesReason] = [
            invalid_ranges_reason
            for invalid_ranges_reason_class in sorted(InvalidRangesReason.__subclasses__(), key=lambda x: x.__name__)
            if (invalid_ranges_reason := invalid_ranges_reason_class.check_and_possibly_return_invalid_part(df))
        ]
        return invalid_ranges_reasons or None


class StartsOrEndsMissingValues(InvalidRangesReason):
    @property
    def reason(self) -> str:  # noqa: D102
        return f"{self.invalid_part.shape[0]} starts or ends are nan"

    @staticmethod
    def check_and_possibly_return_invalid_part(df: "pd.DataFrame") -> "InvalidRangesReason | None":  # noqa: D102
        any_missing_values = df[RANGE_COLS].isna()
        if any_missing_values.any().any():
            return StartsOrEndsMissingValues(df.loc[any_missing_values.any(axis=1)])
        return None


class EmptyOrNegativeIntervals(InvalidRangesReason):
    @property
    def reason(self) -> str:  # noqa: D102
        return f"{self.invalid_part.shape[0]} intervals are empty or negative length (end <= start)"

    @staticmethod
    def check_and_possibly_return_invalid_part(df: "pd.DataFrame") -> "InvalidRangesReason | None":  # noqa: D102
        if np.any(empty_or_negative_filter := (df[START_COL] >= df[END_COL])):
            return EmptyOrNegativeIntervals(df.loc[empty_or_negative_filter])
        return None


class StartOrEndsNegative(InvalidRangesReason):
    @property
    def reason(self) -> str:  # noqa: D102
        return f"{self.invalid_part.shape[0]} starts or ends are < 0"

    @staticmethod
    def check_and_possibly_return_invalid_part(df: "pd.DataFrame") -> "InvalidRangesReason | None":  # noqa: D102
        if (negative_filter := (df[START_COL] < 0) | (df[END_COL] < 0)).any():
            return StartOrEndsNegative(df.loc[negative_filter])
        return None
