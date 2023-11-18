import shutil
from functools import cached_property
from typing import Literal, Iterable

import pandas as pd
from tabulate import tabulate

from pyranges.methods.intersection import _overlap
from pyranges.names import RANGE_COLS
from pyranges.tostring import adjust_table_width


class RangeFrame(pd.DataFrame):
    """Class for range based operations"""
    @cached_property
    def _required_columns(self) -> Iterable[str]:
        return RANGE_COLS[:]

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        missing_columns = set(self._required_columns) - set(self.columns)
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

    def __str__(self, max_col_width: int | None = None, max_total_width: int | None = None) -> str:
        """Return string representation."""

        if len(self) >= 8:
            head = [list(v) for _, v in self.head(4).iterrows()]
            tail = [list(v) for _, v in self.tail(4).iterrows()]
            data = head + [["..."] * self.shape[1]] + tail if len(self) > 8 else head + tail
        else:
            data = [list(v) for _, v in self.iterrows()]

        adjusted_data = adjust_table_width(
            data=data,
            headers=list(self.columns),
            dtypes=[str(t) for t in self.dtypes],
            max_col_width=max_col_width,
            max_total_width=shutil.get_terminal_size().columns if max_total_width is None else max_total_width,
        )
        columns_not_shown = "."
        truncated_data = adjusted_data.truncated_data
        truncated_headers = adjusted_data.truncated_headers
        truncated_dtypes = adjusted_data.truncated_dtypes
        if not len(adjusted_data.truncated_headers) == len(self.columns):
            num_not_shown = len(self.columns) - len(truncated_headers)
            not_shown = [f'"{e}"' for e in self.columns[adjusted_data.included_columns:adjusted_data.included_columns + 3]]
            if num_not_shown > 3:
                not_shown.append("...")
            columns_not_shown = f" ({num_not_shown} columns not shown: {", ".join(not_shown)})."
            truncated_data = [row + ["..."] for row in truncated_data]
            truncated_headers += ["..."]
            truncated_dtypes += ["..."]
        headers_with_dtype = [f"{h}\n{d}" for h, d in zip(truncated_headers, truncated_dtypes)]
        class_and_shape_info = f"{self.__class__.__name__} with {self.shape[0]} rows and {self.shape[1]} columns"
        return f"{tabulate(truncated_data, headers_with_dtype)}\n{class_and_shape_info}{columns_not_shown}"

    def __repr__(self, max_col_width: int | None = None, max_total_width: int | None = None) -> str:
        return self.__str__(max_col_width=max_col_width, max_total_width=max_total_width)

    def overlap(
        self,
        other: "RangeFrame",
        how: Literal["first", "containment", "all"] = "first",
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

        Returns
        -------
        PyRanges
            PyRanges with overlapping ranges.

        Examples
        --------
        >>> import pyranges as pr
        >>> r = pr.RangeFrame({"Start": [1, 1, 2, 2], "End": [3, 3, 5, 4], "Id": list("abcd")})
        >>> r2 = pr.RangeFrame({"Start": [0, 2], "End": [1, 20]})
        >>> r.overlap(r2, how="first")
          Start      End  Id
          int64    int64  object
        -------  -------  --------
              1        3  a
              1        3  b
              2        5  c
              2        4  d

        >>> r.overlap(r2, how="containment")
          Start      End  Id
          int64    int64  object
        -------  -------  --------
              2        5  c
              2        4  d

        >>> r.overlap(r2, how="all")
          Start      End  Id
          int64    int64  object
        -------  -------  --------
              1        3  a
              1        3  b
              2        5  c
              2        4  d
        """
        return RangeFrame(_overlap(self, other, how))
