import inspect
from collections.abc import Callable, Iterable
from typing import Any, TypeVar

import numpy as np
import pandas as pd
import ruranges

from pyranges.core.names import (
    BY_ENTRY_IN_KWARGS,
    END_COL,
    JOIN_SUFFIX,
    RANGE_COLS,
    SKIP_IF_DF_EMPTY_TYPE,
    SKIP_IF_EMPTY_ANY,
    SKIP_IF_EMPTY_BOTH,
    SKIP_IF_EMPTY_LEFT,
    SKIP_IF_EMPTY_RIGHT,
    START_COL,
    VALID_BY_TYPES,
    VALID_COMBINE_OPTIONS,
    VALID_DIRECTION_TYPE,
    VALID_JOIN_TYPE,
    VALID_OVERLAP_TYPE,
    CombineIntervalColumnsOperation,
)
from pyranges.core.pyranges_helpers import arg_to_list, factorize, factorize_binary
from pyranges.core.tostring import tostring
from pyranges.methods.complement_overlaps import _complement_overlaps
from pyranges.methods.join import _both_dfs
from pyranges.methods.merge import _merge
from pyranges.methods.sort import sort_factorize_dict
from pyranges.range_frame.range_frame_validator import InvalidRangesReason


def should_skip_operation(df: pd.DataFrame, *, df2: pd.DataFrame, skip_if_empty: SKIP_IF_DF_EMPTY_TYPE) -> bool:
    """Whether to skip operation because one or more dfs are empty."""
    if df.empty and df2.empty:
        return skip_if_empty in {SKIP_IF_EMPTY_BOTH, SKIP_IF_EMPTY_ANY, SKIP_IF_EMPTY_LEFT, SKIP_IF_EMPTY_RIGHT}
    if df.empty:
        return skip_if_empty in {SKIP_IF_EMPTY_LEFT, SKIP_IF_EMPTY_ANY}
    if df2.empty:
        return skip_if_empty in {SKIP_IF_EMPTY_RIGHT, SKIP_IF_EMPTY_ANY}
    return False


TRangeFrame = TypeVar("TRangeFrame", bound="RangeFrame")


class RangeFrame(pd.DataFrame):
    """Class for range based operations.

    A table with Start and End columns. Parent class of PyRanges. Subclass of pandas DataFrame.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

    def __str__(
        self,
        **kwargs: int | None,
    ) -> str:  # , max_col_width: int | None = None, max_total_width: int | None = None) -> str:
        """Return string representation."""
        str_repr = tostring(
            self,
            max_col_width=kwargs.get("max_col_width"),
            max_total_width=kwargs.get("max_total_width"),
        )
        if reasons := InvalidRangesReason.formatted_reasons_list(self):
            return f"{str_repr}\nInvalid ranges:\n{reasons}"
        return str_repr

    def __repr__(self, max_col_width: int | None = None, max_total_width: int | None = None) -> str:
        return self.__str__(max_col_width=max_col_width, max_total_width=max_total_width)

    def merge_overlaps(
        self,
        *,
        count_col: str | None = None,
        match_by: VALID_BY_TYPES = None,
        slack: int = 0,
    ) -> "RangeFrame":
        match_by = arg_to_list(match_by)
        return _merge(self, by=match_by, count_col=count_col, slack=slack)

    def count_overlaps(
        self,
        other: "RangeFrame",
        *,
        match_by: str | list[str] | None = None,
        slack: int = 0,
    ) -> "pd.Series":
        f1, f2 = factorize_binary(self, other, match_by)
        import ruranges

        return ruranges.count_overlaps_numpy(
            f1,
            self[START_COL].to_numpy(),
            self[END_COL].to_numpy(),
            f2,
            other[START_COL].to_numpy(),
            other[END_COL].to_numpy(),
            slack=slack,
        )

    def combine_interval_columns(
        self,
        function: VALID_COMBINE_OPTIONS | CombineIntervalColumnsOperation = "intersect",
        *,
        start: str = START_COL,
        end: str = END_COL,
        start2: str = START_COL + JOIN_SUFFIX,
        end2: str = END_COL + JOIN_SUFFIX,
        drop_old_columns: bool = True,
    ) -> "RangeFrame":
        """Use two pairs of columns representing intervals to create a new start and end column.

        The function is designed as post-processing after join_ranges to aggregate the coordinates of the two intervals.
        By default, the new start and end columns will be the intersection of the intervals.

        Parameters
        ----------
        function : {"intersect", "union", "swap"} or Callable, default "intersect"
            How to combine the self and other intervals: "intersect", "union", or "swap"
            If a callable is passed, it should take four Series arguments: start1, end1, start2, end2;
            and return a tuple of two integers: (new_starts, new_ends).

        start : str, default "Start"
            Column name for Start of first interval
        end : str, default "End"
            Column name for End of first interval
        start2 : str, default "Start_b"
            Column name for Start of second interval
        end2 : str, default "End_b"
            Column name for End of second interval
        drop_old_columns : bool, default True
            Whether to drop the above mentioned columns.

        Examples
        --------
        >>> gr1, gr2 = pr.example_data.aorta.head(3).remove_nonloc_columns(), pr.example_data.aorta2.head(3).remove_nonloc_columns()
        >>> j = gr1.join_ranges(gr2)
        >>> j
          index  |    Chromosome      Start      End  Strand        Start_b    End_b
          int64  |    category        int64    int64  category        int64    int64
        -------  ---  ------------  -------  -------  ----------  ---------  -------
              0  |    chr1             9916    10115  -                9988    10187
              1  |    chr1             9916    10115  -               10079    10278
              2  |    chr1             9939    10138  +               10073    10272
              3  |    chr1             9951    10150  -                9988    10187
              4  |    chr1             9951    10150  -               10079    10278
        PyRanges with 5 rows, 6 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        The default operation is to intersect the intervals:

        >>> j.combine_interval_columns()
          index  |    Chromosome      Start      End  Strand
          int64  |    category        int64    int64  category
        -------  ---  ------------  -------  -------  ----------
              0  |    chr1             9988    10115  -
              1  |    chr1            10079    10115  -
              2  |    chr1            10073    10138  +
              3  |    chr1             9988    10150  -
              4  |    chr1            10079    10150  -
        PyRanges with 5 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        Take the union instead:

        >>> j.combine_interval_columns('union')
          index  |    Chromosome      Start      End  Strand
          int64  |    category        int64    int64  category
        -------  ---  ------------  -------  -------  ----------
              0  |    chr1             9916    10187  -
              1  |    chr1             9916    10278  -
              2  |    chr1             9939    10272  +
              3  |    chr1             9951    10187  -
              4  |    chr1             9951    10278  -
        PyRanges with 5 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        >>> j.combine_interval_columns('swap')
          index  |    Chromosome      Start      End  Strand
          int64  |    category        int64    int64  category
        -------  ---  ------------  -------  -------  ----------
              0  |    chr1             9988    10187  -
              1  |    chr1            10079    10278  -
              2  |    chr1            10073    10272  +
              3  |    chr1             9988    10187  -
              4  |    chr1            10079    10278  -
        PyRanges with 5 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.


        Use a custom function that keeps the start of the first interval and the end of the second:

        >>> def custom_combine(s1, e1, s2, e2): return (s1, e2)
        >>> j.combine_interval_columns(custom_combine)
          index  |    Chromosome      Start      End  Strand
          int64  |    category        int64    int64  category
        -------  ---  ------------  -------  -------  ----------
              0  |    chr1             9916    10187  -
              1  |    chr1             9916    10278  -
              2  |    chr1             9939    10272  +
              3  |    chr1             9951    10187  -
              4  |    chr1             9951    10278  -
        PyRanges with 5 rows, 4 columns, and 1 index columns.
        Contains 1 chromosomes and 2 strands.

        """
        from pyranges.methods.combine_positions import (
            _intersect_interval_columns,
            _swap_interval_columns,
            _union_interval_columns,
        )

        if function == "intersect":
            function = _intersect_interval_columns
        elif function == "union":
            function = _union_interval_columns
        elif function == "swap":
            function = _swap_interval_columns

        new_starts, new_ends = function(self[start], self[end], self[start2], self[end2])

        z = self.copy()
        z[START_COL] = new_starts
        z[END_COL] = new_ends

        cols_to_drop = list({start, end, start2, end2}.difference(RANGE_COLS) if drop_old_columns else {})

        return z.drop_and_return(labels=cols_to_drop, axis="columns")

    def cluster(
        self,
        *,
        match_by: VALID_BY_TYPES = None,
        cluster_column: str = "Cluster",
        slack: int = 0,
    ) -> "RangeFrame":
        match_by = arg_to_list(match_by)

        factorized = factorize(self, match_by)
        cluster, idx = ruranges.cluster_numpy(
            factorized,
            self[START_COL].to_numpy(),
            self[END_COL].to_numpy(),
            slack,
        )

        res = self.take(idx).copy()
        res.insert(res.shape[1], cluster_column, cluster)
        return _mypy_ensure_rangeframe(res)

    def complement_overlaps(
        self: "RangeFrame",
        other: "RangeFrame",
        *,
        match_by: VALID_BY_TYPES = None,
        slack: int = 0,
    ) -> "RangeFrame":
        match_by = arg_to_list(match_by)
        return _complement_overlaps(self, other, by=match_by, slack=slack)

    def join_ranges(
        self,
        other: "RangeFrame",
        *,
        join_type: VALID_JOIN_TYPE = "inner",
        multiple: VALID_OVERLAP_TYPE = "all",
        match_by: VALID_BY_TYPES = None,
        slack: int = 0,
        suffix: str = JOIN_SUFFIX,
        contained_intervals_only: bool = False,
        report_overlap_column: str | None = None,
    ) -> "RangeFrame":
        res = _both_dfs(
            self,
            other,
            by=match_by,
            slack=slack,
            multiple=multiple,
            contained=contained_intervals_only,
            join_type=join_type,
            suffix=suffix,
        )

        if report_overlap_column:
            res[report_overlap_column] = res[["End", "End" + suffix]].min(axis=1) - res[
                ["Start", "Start" + suffix]
            ].max(axis=1)

        res.index = res.index.astype(np.int64)

        return _mypy_ensure_rangeframe(res)

    def max_disjoint(
        self,
        *,
        slack: int = 0,
        match_by: VALID_BY_TYPES = None,
    ) -> "RangeFrame":
        import ruranges

        factorized = factorize(self, match_by)

        idx = ruranges.max_disjoint_numpy(
            factorized,
            self[START_COL].to_numpy(),
            self[END_COL].to_numpy(),
            slack,
        )
        return _mypy_ensure_rangeframe(self.take(idx))

    def nearest(
        self,
        other: "RangeFrame",
        *,
        match_by: VALID_BY_TYPES = None,
        suffix: str = JOIN_SUFFIX,
        exclude_overlaps: bool = False,
        k: int = 1,
        dist_col: str | None = "Distance",
        direction: VALID_DIRECTION_TYPE = "any",
    ) -> "RangeFrame":
        f1, f2 = factorize_binary(self, other, match_by)
        idx1, idx2, dist = ruranges.nearest_numpy(
            chrs=f1.astype(np.uint32),
            starts=self[START_COL].to_numpy(),
            ends=self[END_COL].to_numpy(),
            chrs2=f2.astype(np.uint32),
            starts2=other[START_COL].to_numpy(),
            ends2=other[END_COL].to_numpy(),
            k=k,
            slack=0,
            include_overlaps=not exclude_overlaps,
            direction=direction,
        )

        left = self.take(idx1)
        to_concat = [
            left.reset_index(drop=True),
            pd.DataFrame(other).take(idx2).add_suffix(suffix).reset_index(drop=True),
        ]
        if dist_col is not None:
            to_concat.append(pd.Series(dist, name=dist_col))

        res = pd.concat(to_concat, axis=1)
        res.index = pd.Index(left.index.to_numpy())

        return _mypy_ensure_rangeframe(res)

    def overlap(
        self,
        other: "RangeFrame",
        multiple: VALID_OVERLAP_TYPE = "all",
        slack: int = 0,
        *,
        contained_intervals_only: bool = False,
        match_by: VALID_BY_TYPES = None,
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

        >>> r.overlap(r2, multiple="first")
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              0  |          1        3  a
              1  |          1        3  b
              3  |          2        4  d
              2  |          2        5  a
        RangeFrame with 4 rows, 3 columns, and 1 index columns.

        >>> r.overlap(r2, contained_intervals_only=True)
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              3  |          2        4  d
              2  |          2        5  a
        RangeFrame with 2 rows, 3 columns, and 1 index columns.

        >>> r.overlap(r2, multiple="all")
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              0  |          1        3  a
              1  |          1        3  b
              3  |          2        4  d
              2  |          2        5  a
        RangeFrame with 4 rows, 3 columns, and 1 index columns.

        >>> r.overlap(r2, multiple="all", match_by="Id")
          index  |      Start      End  Id
          int64  |      int64    int64  object
        -------  ---  -------  -------  --------
              3  |          2        4  d
        RangeFrame with 1 rows, 3 columns, and 1 index columns.

        """
        from pyranges.methods.overlap import _overlap

        by = arg_to_list(match_by)

        result = _overlap(
            self,
            other,
            by=by,
            slack=slack,
            multiple=multiple,
            contained=contained_intervals_only,
        )

        return _mypy_ensure_rangeframe(result)

    def sort_ranges(
        self: "RangeFrame",
        match_by: VALID_BY_TYPES = None,
        *,
        natsort: bool = True,
    ):
        by = arg_to_list(match_by)
        by_sort_order_as_int = sort_factorize_dict(self, by, use_natsort=natsort)
        idxs = ruranges.sort_intervals_numpy(
            by_sort_order_as_int,
            self[START_COL].to_numpy(),
            self[END_COL].to_numpy(),
        )
        return self.take(idxs)

    def subtract_ranges(
        self: "RangeFrame",
        other: "RangeFrame",
        match_by: VALID_BY_TYPES = None,
    ) -> "RangeFrame":
        f1, f2 = factorize_binary(self, other, match_by)

        idx, start, end = ruranges.subtract_numpy(
            f1,
            self[START_COL].to_numpy(),
            self[END_COL].to_numpy(),
            f2,
            other[START_COL].to_numpy(),
            other[END_COL].to_numpy(),
        )

        output = self.take(idx).copy()
        output[START_COL], output[END_COL] = start, end
        return _mypy_ensure_rangeframe(output)

    def sort_by_position(self) -> "RangeFrame":
        """Sort by Start and End columns."""
        return _mypy_ensure_rangeframe(self.sort_values(RANGE_COLS))

    def reasons_why_frame_is_invalid(self) -> list[InvalidRangesReason] | None:  # noqa: D102
        __doc__ = InvalidRangesReason.is_invalid_ranges_reasons.__doc__  # noqa: A001, F841

        return InvalidRangesReason.is_invalid_ranges_reasons(self)

    def copy(self, *args, **kwargs) -> "RangeFrame":  # noqa: D102
        return _mypy_ensure_rangeframe(super().copy(*args, **kwargs))

    def drop(self, *args, **kwargs) -> "RangeFrame | None":  # type: ignore[override]  # noqa: D102
        return self.__class__(super().drop(*args, **kwargs))

    def drop_and_return[T: "RangeFrame"](self: T, *args: Any, **kwargs: Any) -> T:  # noqa: PYI019, D102
        kwargs["inplace"] = False
        return self.__class__(super().drop(*args, **kwargs))

    def reindex(self, *args, **kwargs) -> "RangeFrame":  # noqa: D102
        return self.__class__(super().reindex(*args, **kwargs))


def _mypy_ensure_rangeframe(r: pd.DataFrame) -> "RangeFrame":
    result = RangeFrame(r)
    if not isinstance(result, RangeFrame):
        msg = f"Expected RangeFrame, got {type(result)}"
        raise TypeError(msg)
    return result


def assert_valid_ranges(function: Callable, *args: "RangeFrame") -> None:
    """Raise ValueError because function is called on invalid ranges."""
    if any(r.reasons_why_frame_is_invalid() for r in args):
        is_not_lambda = function.__name__ != "<lambda>"
        function_repr = function.__name__ if is_not_lambda else inspect.getsource(function).strip()
        msg = f"Cannot perform function on invalid ranges (function was {function_repr})."
        raise ValueError(msg)


def _with_group_keys_to_kwargs(by: list[str]) -> Callable:
    def decorator(func) -> Callable:
        def wrapper(group, **kwargs) -> Any:
            # Extract the group keys and add them to kwargs
            names = group.name if isinstance(group.name, Iterable) and not isinstance(group.name, str) else [group.name]
            kwargs[BY_ENTRY_IN_KWARGS] = dict(zip(by, names, strict=True))
            return func(group, **kwargs)

        return wrapper

    return decorator
