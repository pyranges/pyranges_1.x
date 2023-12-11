from collections.abc import Callable, Iterable
from functools import cached_property
from typing import Any

import pandas as pd

from pyranges.methods.overlap import _overlap
from pyranges.names import RANGE_COLS, VALID_BY_TYPES, VALID_OVERLAP_TYPE
from pyranges.tostring import tostring


class ColUpdater:
    r: "RangeFrame"

    def __init__(self, r: "RangeFrame") -> None:
        self.__dict__["r"] = r

    def __setattr__(
        self,
        key: str,
        value: Any,
    ) -> None:
        if key not in self.r:
            self.r.insert(self.r.shape[-1], key, value)
        else:
            self.r.loc[:, key] = value

    def __setitem__(self, key: str, value: Any) -> None:
        if isinstance(value, pd.Series):
            value = value.to_numpy()
        if key not in self.r.columns:
            self.r.insert(self.r.shape[-1], key, value)
        else:
            self.r.loc[:, key] = value


class RangeFrame(pd.DataFrame):
    """Class for range based operations."""

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

    @property
    def col(self) -> "ColUpdater":
        """Add or update a column.

        Returns
        -------
        ColUpdater

        Examples
        --------
        >>> gr = RangeFrame({"Start": [0, 1], "End": [2, 3]})
        >>> gr.col.Frame = ["Hi", "There"]
        >>> gr
          index  |      Start      End  Frame
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              0  |          0        2  Hi
              1  |          1        3  There
        DataFrame with 2 rows, 3 columns, and 1 index columns.

        >>> gr.col.Frame = ["Next", "words"]
        >>> gr
          index  |      Start      End  Frame
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              0  |          0        2  Next
              1  |          1        3  words
        DataFrame with 2 rows, 3 columns, and 1 index columns.

        >>> gr.col["Start"] = [5000, 1000]
        >>> gr
          index  |      Start      End  Frame
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              0  |       5000        2  Next
              1  |       1000        3  words
        DataFrame with 2 rows, 3 columns, and 1 index columns.
        """
        return ColUpdater(self)

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
        PyRanges
            PyRanges with overlapping ranges.

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
        DataFrame with 4 rows, 3 columns, and 1 index columns.

        >>> r2 = pr.RangeFrame({"Start": [0, 2], "End": [1, 20], "Id": list("ad")})
        >>> r2
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              0  |          0        1  a
              1  |          2       20  d
        DataFrame with 2 rows, 3 columns, and 1 index columns.

        >>> r.overlap(r2, how="first")
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              0  |          1        3  a
              1  |          1        3  b
              2  |          2        5  a
              3  |          2        4  d
        DataFrame with 4 rows, 3 columns, and 1 index columns.

        >>> r.overlap(r2, how="containment")
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              2  |          2        5  a
              3  |          2        4  d
        DataFrame with 2 rows, 3 columns, and 1 index columns.

        >>> r.overlap(r2, how="all")
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              0  |          1        3  a
              1  |          1        3  b
              2  |          2        5  a
              3  |          2        4  d
        DataFrame with 4 rows, 3 columns, and 1 index columns.

        >>> r.overlap(r2, how="all", by="Id")
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              3  |          2        4  d
        DataFrame with 1 rows, 3 columns, and 1 index columns.
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


def _mypy_ensure_rangeframe(r: pd.DataFrame) -> "RangeFrame":
    return RangeFrame(r)
