import abc
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from pyranges.names import END_COL, RANGE_COLS, START_COL

if TYPE_CHECKING:
    import pandas as pd


class InvalidRangesReason:
    """Describe why a range is invalid."""

    @property
    @abc.abstractmethod
    def invalid_part(self) -> "pd.DataFrame":
        return self._invalid_part

    @invalid_part.setter
    def invalid_part(self, df: "pd.DataFrame") -> None:
        self._invalid_part = df

    @property
    @abc.abstractmethod
    def reason(self) -> str:
        pass

    @staticmethod
    @abc.abstractmethod
    def check_and_possibly_return_invalid_part(df: "pd.DataFrame") -> "InvalidRangesReason | None":
        pass

    def to_dict(self) -> dict[str, "str | pd.DataFrame"]:
        return {"reason": self.reason, "invalid_part": self.invalid_part}

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
        >>> invalid = RangeFrame({"Start": [-5, 1, 100, np.nan], "End": [3, 2, 10, 40]})
        >>> invalid
          index  |        Start      End
          int64  |      float64    int64
        -------  ---  ---------  -------
              0  |           -5        3
              1  |            1        2
              2  |          100       10
              3  |          nan       40
        RangeFrame with 4 rows, 2 columns, and 1 index columns.
        >>> invalid_ranges_reasons = invalid.reasons_why_frame_is_invalid()
        >>> for invalid_range in invalid_ranges_reasons:
        ...     print(invalid_range.reason)
        ...     print(invalid_range.invalid_part)
        ...     print()
        Some starts or ends are nan.
          index  |        Start      End
          int64  |      float64    int64
        -------  ---  ---------  -------
              3  |          nan       40
        RangeFrame with 1 rows, 2 columns, and 1 index columns.
        <BLANKLINE>
        Some intervals are empty or negative length (end <= start).
          index  |        Start      End
          int64  |      float64    int64
        -------  ---  ---------  -------
              2  |          100       10
        RangeFrame with 1 rows, 2 columns, and 1 index columns.
        <BLANKLINE>
        Some starts or ends are < 0.
          index  |        Start      End
          int64  |      float64    int64
        -------  ---  ---------  -------
              0  |           -5        3
        RangeFrame with 1 rows, 2 columns, and 1 index columns.
        <BLANKLINE>
        """
        invalid_ranges_reasons: list["InvalidRangesReason"] = [
            invalid_ranges_reason
            for invalid_ranges_reason_class in InvalidRangesReason.__subclasses__()
            if (invalid_ranges_reason := invalid_ranges_reason_class.check_and_possibly_return_invalid_part(df))
        ]
        return invalid_ranges_reasons or None


@dataclass
class StartsOrEndsMissingValues(InvalidRangesReason):
    invalid_part: "pd.DataFrame"
    reason: str = "Some starts or ends are nan."

    @staticmethod
    def check_and_possibly_return_invalid_part(df: "pd.DataFrame") -> "InvalidRangesReason | None":
        any_missing_values = df[RANGE_COLS].isna()
        if any_missing_values.any().any():
            return StartsOrEndsMissingValues(df.loc[any_missing_values.any(axis=1)])
        return None


@dataclass
class EmptyOrNegativeIntervals(InvalidRangesReason):
    invalid_part: "pd.DataFrame"
    reason: str = "Some intervals are empty or negative length (end <= start)."

    @staticmethod
    def check_and_possibly_return_invalid_part(df: "pd.DataFrame") -> "InvalidRangesReason | None":
        if np.any(empty_or_negative_filter := (df[START_COL] >= df[END_COL])):
            return EmptyOrNegativeIntervals(df.loc[empty_or_negative_filter])
        return None


@dataclass
class StartOrEndsNegative(InvalidRangesReason):
    invalid_part: "pd.DataFrame"
    reason: str = "Some starts or ends are < 0."

    @staticmethod
    def check_and_possibly_return_invalid_part(df: "pd.DataFrame") -> "InvalidRangesReason | None":
        if np.any(negative_filter := (df[RANGE_COLS] < 0)):
            return StartOrEndsNegative(df.loc[negative_filter.any(axis="columns")])
        return None
