from collections.abc import Callable, Iterable
from functools import cached_property

import pandas as pd

from pyranges.methods.overlap import _overlap
from pyranges.names import RANGE_COLS, VALID_BY_TYPES, VALID_OVERLAP_TYPE
from pyranges.range_frame.range_frame_validator import InvalidRangesReason
from pyranges.tostring import tostring


class RangeFrame(pd.DataFrame):
    """Class for range based operations."""

    def __new__(cls, *args, **kwargs) -> "RangeFrame | pd.DataFrame":  # type: ignore[misc]
        """Create a new instance of a PyRanges object."""
        # __new__ is a special static method used for creating and
        # returning a new instance of a class. It is called before
        # __init__ and is typically used in scenarios requiring
        # control over the creation of new instances

        if not args:
            return super().__new__(cls)
        if not kwargs:
            return super().__new__(cls)

        df = pd.DataFrame(kwargs.get("data") or (args[0]))

        missing_any_required_columns = not set(RANGE_COLS).issubset(df.columns)
        if missing_any_required_columns:
            return df

        return super().__new__(cls)

    @property
    def _constructor(self) -> type:
        return RangeFrame

    @cached_property
    def _required_columns(self) -> Iterable[str]:
        return RANGE_COLS[:]

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

        missing_columns = [c for c in RANGE_COLS if c not in self.columns]
        if missing_columns:
            msg = f"Missing required columns: {missing_columns}"
            raise ValueError(msg)

    def __str__(
        self,
        **kwargs,
    ) -> str:  # , max_col_width: int | None = None, max_total_width: int | None = None) -> str:
        """Return string representation."""
        return tostring(self, max_col_width=kwargs.get("max_col_width"), max_total_width=kwargs.get("max_total_width"))

    def __repr__(self, max_col_width: int | None = None, max_total_width: int | None = None) -> str:
        return self.__str__(max_col_width=max_col_width, max_total_width=max_total_width)

    def overlap(
        self,
        other: "RangeFrame",
        how: VALID_OVERLAP_TYPE = "first",
        by: VALID_BY_TYPES = None,
        **_,
    ) -> "RangeFrame":
        """Find intervals in self overlapping other..

        Parameters
        ----------
        other
            Other ranges to find overlaps with.
        how
            How to find overlaps. "first" finds the first overlap, "containment" finds all overlaps
            where self is contained in other, and "all" finds all overlaps.
        by:
            Grouping columns. If None, all columns are used.

        Returns
        -------
        RangeFrame
            RangeFrame with overlapping ranges.

        Examples
        --------
        >>> import pyranges as pr
        >>> r = pr.RangeFrame({"Start": [1, 1, 2, 2], "End": [3, 3, 5, 4], "Id": list("abad")})
        >>> r
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              0  |          1        3  a
              1  |          1        3  b
              2  |          2        5  a
              3  |          2        4  d
        RangeFrame with 4 rows, 3 columns, and 1 index columns.

        >>> r2 = pr.RangeFrame({"Start": [0, 2], "End": [1, 20], "Id": list("ad")})
        >>> r2
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              0  |          0        1  a
              1  |          2       20  d
        RangeFrame with 2 rows, 3 columns, and 1 index columns.

        >>> r.overlap(r2, how="first")
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              0  |          1        3  a
              1  |          1        3  b
              2  |          2        5  a
              3  |          2        4  d
        RangeFrame with 4 rows, 3 columns, and 1 index columns.

        >>> r.overlap(r2, how="containment")
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              2  |          2        5  a
              3  |          2        4  d
        RangeFrame with 2 rows, 3 columns, and 1 index columns.

        >>> r.overlap(r2, how="all")
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              0  |          1        3  a
              1  |          1        3  b
              2  |          2        5  a
              3  |          2        4  d
        RangeFrame with 4 rows, 3 columns, and 1 index columns.

        >>> r.overlap(r2, how="all", by="Id")
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              3  |          2        4  d
        RangeFrame with 1 rows, 3 columns, and 1 index columns.
        """
        return self.apply_pair(other, _overlap, how=how, by=by)

    def apply_single(
        self,
        function: Callable[["RangeFrame"], "RangeFrame"],
        by: VALID_BY_TYPES,
        **kwargs,
    ) -> "RangeFrame":
        """Call a function on a RangeFrame.

        Parameters
        ----------
        function: Callable
            Function to call.

        by: str or list of str
            Group by these columns.

        kwargs: dict
            Passed to function.
        """
        if not by:
            return RangeFrame(function(self, **kwargs))
        by = self._by_to_list(by)
        return _mypy_ensure_rangeframe(
            self.groupby(by).apply(function, by=by, **kwargs).reset_index(drop=True),  # type: ignore[arg-type]
        )

    def apply_pair(
        self,
        other: "RangeFrame",
        function: Callable[["RangeFrame", "RangeFrame"], "RangeFrame"],
        by: VALID_BY_TYPES = None,
        **kwargs,
    ) -> "RangeFrame":
        """Call a function on two RangeFrames.

        Parameters
        ----------
        other: RangeFrame
            Other RangeFrame.

        function: Callable
            Function to call.

        by: str or list of str, default None
            Group by these columns.

        kwargs: dict
            Passed to function.
        """
        if by is None:
            return RangeFrame(function(self, other, **kwargs))

        by = self._by_to_list(by)
        results = []
        empty = RangeFrame(columns=other.columns)
        others = dict(list(other.groupby(by)))

        for key, _df in self.groupby(by):
            odf = others.get(key, empty)
            results.append(function(_mypy_ensure_rangeframe(_df), _mypy_ensure_rangeframe(odf), **kwargs))

        return RangeFrame(pd.concat(results))

    def sort_by_position(self) -> "RangeFrame":
        """Sort by Start and End columns."""
        return _mypy_ensure_rangeframe(self.sort_values(RANGE_COLS))

    @staticmethod
    def _by_to_list(by: str | Iterable[str] | None) -> list[str]:
        return [by] if isinstance(by, str) else ([*by] if by is not None else [])

    def reasons_why_frame_is_invalid(self) -> list[InvalidRangesReason] | None:
        __doc__ = InvalidRangesReason.is_invalid_ranges_reasons.__doc__  # noqa: A001, F841

        return InvalidRangesReason.is_invalid_ranges_reasons(self)


def _mypy_ensure_rangeframe(r: pd.DataFrame) -> "RangeFrame":
    return RangeFrame(r)
