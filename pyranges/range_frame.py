from dataclasses import dataclass
from functools import cached_property
from typing import Literal, Iterable, Self

import pandas as pd

from pyranges.methods.overlap import _overlap
from pyranges.names import RANGE_COLS
from pyranges.tostring import tostring


class ColUpdater:
    def __init__(self, r: "RangeFrame") -> Self:
        self.__dict__["r"] = r

    def __setattr__(self, key, value) -> None:
        if key not in self.r:
            self.r.insert(self.r.shape[-1], key, value)
        else:
            self.r.loc[:, key] = value

    def __setitem__(self, key: str, value: Iterable) -> None:
        if key not in self.r:
            self.r.insert(self.r.shape[-1], key, value)
        else:
            self.r.loc[:, key] = value


class RangeFrame(pd.DataFrame):
    """Class for range based operations"""

    @cached_property
    def _required_columns(self) -> Iterable[str]:
        return RANGE_COLS[:]

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        missing_columns = [c for c in self._required_columns if c not in self.columns]
        if missing_columns:
            msg = f"Missing required columns: {missing_columns}"
            raise ValueError(msg)
        self._check_index_column_names()

    def _check_index_column_names(self):
        # Check for columns with the same name as the index
        if self.index.name is not None and self.index.name in self.columns:
            raise ValueError(
                f"Index name '{self.index.name}' cannot be the same as a column name."
            )
        if isinstance(self.index, pd.MultiIndex):
            for level in self.index.names:
                if level in self.columns:
                    raise ValueError(
                        f"Index level name '{level}' cannot be the same as a column name."
                    )

    def __str__(
        self, max_col_width: int | None = None, max_total_width: int | None = None
    ) -> str:
        """Return string representation."""
        return tostring(
            self, max_col_width=max_col_width, max_total_width=max_total_width
        )

    def __repr__(
        self, max_col_width: int | None = None, max_total_width: int | None = None
    ) -> str:
        return self.__str__(
            max_col_width=max_col_width, max_total_width=max_total_width
        )

    @property
    def col(self) -> "ColUpdater":
        """Add or update a column.

        Returns
        -------

        Examples
        >>> gr = RangeFrame({"Start": [0, 1], "End": [2, 3]})
        >>> gr.col.Frame = ["Hi", "There"]
        >>> gr
          Start      End  Frame
          int64    int64  object
        -------  -------  --------
              0        2  Hi
              1        3  There
        RangeFrame with 2 rows and 3 columns.

        >>> gr.col.Frame = ["Next", "words"]
        >>> gr
          Start      End  Frame
          int64    int64  object
        -------  -------  --------
              0        2  Next
              1        3  words
        RangeFrame with 2 rows and 3 columns.

        >>> gr.col["Start"] = [5000, 1000]
        >>> gr
          Start      End  Frame
          int64    int64  object
        -------  -------  --------
           5000        2  Next
           1000        3  words
        RangeFrame with 2 rows and 3 columns.
        """

        return ColUpdater(self)

    def overlap(
        self,
        other: "RangeFrame",
        how: Literal["first", "containment", "all"] = "first",
        by: str | list[str] | None = None,
        **_,
    ) -> "RangeFrame":
        """
        Find intervals in self overlapping other..

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
          Start      End  Id
          int64    int64  object
        -------  -------  --------
              1        3  a
              1        3  b
              2        5  a
              2        4  d
        RangeFrame with 4 rows and 3 columns.
        >>> r2 = pr.RangeFrame({"Start": [0, 2], "End": [1, 20], "Id": list("ad")})
        >>> r2
          Start      End  Id
          int64    int64  object
        -------  -------  --------
              0        1  a
              2       20  d
        RangeFrame with 2 rows and 3 columns.

        >>> r.overlap(r2, how="first")
          Start      End  Id
          int64    int64  object
        -------  -------  --------
              1        3  a
              1        3  b
              2        5  a
              2        4  d
        RangeFrame with 4 rows and 3 columns.

        >>> r.overlap(r2, how="containment")
          Start      End  Id
          int64    int64  object
        -------  -------  --------
              2        5  a
              2        4  d
        RangeFrame with 2 rows and 3 columns.

        >>> r.overlap(r2, how="all")
          Start      End  Id
          int64    int64  object
        -------  -------  --------
              1        3  a
              1        3  b
              2        5  a
              2        4  d
        RangeFrame with 4 rows and 3 columns.

        >>> r.overlap(r2, how="all", by="Id")
          Start      End  Id
          int64    int64  object
        -------  -------  --------
              2        4  d
        RangeFrame with 1 rows and 3 columns.
        """
        return self.apply_pair(other, _overlap, how=how, by=by)

    def apply_single(
        self, function, by: str | list[str] | None = None, **kwargs,
    ) -> "pr.RangeFrame":
        if by is None:
            return RangeFrame(function(self, **kwargs))
        res = self.groupby(by).apply(function, by=by, **kwargs).reset_index(drop=True)

        return res

    def apply_pair(self, other, function, by, **kwargs):
        if by is None:
            return RangeFrame(function(self, other, **kwargs))

        """
        If opposite

        add temp column to match on

        run super method with added column

        remove temp column afterwards
        """

        results = []
        empty = RangeFrame(columns=other.columns)
        others = {k: v for k, v in other.groupby(by)}

        for key, _df in self.groupby(by):
            odf = others.get(key, empty)
            results.append(function(_df, odf, **kwargs))
        return RangeFrame(pd.concat(results))
