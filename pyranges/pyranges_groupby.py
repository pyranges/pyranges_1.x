from collections.abc import Callable
from typing import TYPE_CHECKING, Any

import pandas as pd
import pandas.core.groupby  # type: ignore[name-defined]

from pyranges.core.names import return_pyranges_if_possible

if TYPE_CHECKING:
    import pyranges as pr


class PyRangesDataFrameGroupBy(pandas.core.groupby.DataFrameGroupBy):
    def __init__(self, pandas_groupby) -> None:
        self._pandas_groupby = pandas_groupby

    @property
    def pandas_groupby(self) -> pd.core.groupby.DataFrameGroupBy:
        """Return the underlying pandas groupby object."""
        return self.__dict__["_pandas_groupby"]

    def __getattr__(self, item) -> Any:
        # Handle attribute access, e.g., g.some_method()
        pd_grpby = self.__dict__["_pandas_groupby"]

        if hasattr(pd_grpby, item):
            attr = getattr(pd_grpby, item)
            if callable(attr):

                def wrapper(*args, **kwargs) -> Callable:
                    result = attr(*args, **kwargs)
                    return return_pyranges_if_possible(result)

                return wrapper
            return attr
        return None

    def __getitem__(self, key) -> "PyRangesDataFrameGroupBy":
        # Handle item access, e.g., g['column_name']
        result = self.pandas_groupby[key]
        return PyRangesDataFrameGroupBy(result)

    @return_pyranges_if_possible
    def agg(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.agg(*args, **kwargs)

    @return_pyranges_if_possible
    def aggregate(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.aggregate(*args, **kwargs)

    @return_pyranges_if_possible
    def all(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.all(*args, **kwargs)

    @return_pyranges_if_possible
    def any(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.any(*args, **kwargs)

    @return_pyranges_if_possible
    def apply(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.apply(*args, **kwargs)

    @return_pyranges_if_possible
    def bfill(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.bfill(*args, **kwargs)

    @return_pyranges_if_possible
    def cumcount(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.cumcount(*args, **kwargs)

    @return_pyranges_if_possible
    def cummax(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.cummax(*args, **kwargs)

    @return_pyranges_if_possible
    def cummin(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.cummin(*args, **kwargs)

    @return_pyranges_if_possible
    def cumprod(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.cumprod(*args, **kwargs)

    @return_pyranges_if_possible
    def cumsum(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.cumsum(*args, **kwargs)

    @return_pyranges_if_possible
    def describe(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.describe(*args, **kwargs)

    @return_pyranges_if_possible
    def diff(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.diff(*args, **kwargs)

    @return_pyranges_if_possible
    def ewm(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.ewm(*args, **kwargs)

    @return_pyranges_if_possible
    def expanding(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.expanding(*args, **kwargs)

    @return_pyranges_if_possible
    def ffill(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.ffill(*args, **kwargs)

    @return_pyranges_if_possible
    def fillna(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.fillna(*args, **kwargs)

    @return_pyranges_if_possible
    def filter(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.filter(*args, **kwargs)

    @return_pyranges_if_possible
    def first(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.first(*args, **kwargs)

    @return_pyranges_if_possible
    def get_group(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.get_group(*args, **kwargs)

    @return_pyranges_if_possible
    def head(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.head(*args, **kwargs)

    @return_pyranges_if_possible
    def idxmax(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.idxmax(*args, **kwargs)

    @return_pyranges_if_possible
    def idxmin(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.idxmin(*args, **kwargs)

    @return_pyranges_if_possible
    def last(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.last(*args, **kwargs)

    @return_pyranges_if_possible
    def max(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.max(*args, **kwargs)

    @return_pyranges_if_possible
    def mean(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.mean(*args, **kwargs)

    @return_pyranges_if_possible
    def median(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.median(*args, **kwargs)

    @return_pyranges_if_possible
    def min(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.min(*args, **kwargs)

    @return_pyranges_if_possible
    def ngroup(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.ngroup(*args, **kwargs)

    @return_pyranges_if_possible
    def nunique(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.nunique(*args, **kwargs)

    @return_pyranges_if_possible
    def ohlc(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.ohlc(*args, **kwargs)

    @return_pyranges_if_possible
    def pct_change(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.pct_change(*args, **kwargs)

    @return_pyranges_if_possible
    def pipe(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.pipe(*args, **kwargs)

    @return_pyranges_if_possible
    def prod(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.prod(*args, **kwargs)

    @return_pyranges_if_possible
    def quantile(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.quantile(*args, **kwargs)

    @return_pyranges_if_possible
    def rank(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.rank(*args, **kwargs)

    @return_pyranges_if_possible
    def resample(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.resample(*args, **kwargs)

    @return_pyranges_if_possible
    def rolling(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.rolling(*args, **kwargs)

    @return_pyranges_if_possible
    def sample(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.sample(*args, **kwargs)

    @return_pyranges_if_possible
    def sem(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.sem(*args, **kwargs)

    @return_pyranges_if_possible
    def shift(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.shift(*args, **kwargs)

    @return_pyranges_if_possible
    def size(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.size(*args, **kwargs)

    @return_pyranges_if_possible
    def skew(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.skew(*args, **kwargs)

    @return_pyranges_if_possible
    def std(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.std(*args, **kwargs)

    @return_pyranges_if_possible
    def sum(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.sum(*args, **kwargs)

    @return_pyranges_if_possible
    def tail(self, *args, **kwargs) -> "pr.PyRanges | pd.DataFrame | pd.Series":  # noqa: D102
        return self.pandas_groupby.tail(*args, **kwargs)
