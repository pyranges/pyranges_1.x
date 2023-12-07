import pandas as pd

import pyranges as pr
from pyranges.names import GENOME_LOC_COLS


def return_pyranges_if_possible(method):
    def wrapper(*args, **kwargs):
        # Call the original groupby method
        result = method(*args, **kwargs)

        # Check if the result should be a MySpecialDataFrame
        if isinstance(result, pd.DataFrame) and set(GENOME_LOC_COLS).issubset(result.columns):
            return pr.PyRanges(result)
        else:
            return result

    return wrapper


class PyRangesGroupBy(pd.core.groupby.DataFrameGroupBy):

    @return_pyranges_if_possible
    def agg(self, *args, **kwargs):
        return super().agg(*args, **kwargs)

    @return_pyranges_if_possible
    def aggregate(self, *args, **kwargs):
        return super().aggregate(*args, **kwargs)

    @return_pyranges_if_possible
    def all(self, *args, **kwargs):
        return super().all(*args, **kwargs)

    @return_pyranges_if_possible
    def any(self, *args, **kwargs):
        return super().any(*args, **kwargs)

    @return_pyranges_if_possible
    def apply(self, *args, **kwargs):
        return super().apply(*args, **kwargs)

    @return_pyranges_if_possible
    def bfill(self, *args, **kwargs):
        return super().bfill(*args, **kwargs)

    @return_pyranges_if_possible
    def cumcount(self, *args, **kwargs):
        return super().cumcount(*args, **kwargs)

    @return_pyranges_if_possible
    def cummax(self, *args, **kwargs):
        return super().cummax(*args, **kwargs)

    @return_pyranges_if_possible
    def cummin(self, *args, **kwargs):
        return super().cummin(*args, **kwargs)

    @return_pyranges_if_possible
    def cumprod(self, *args, **kwargs):
        return super().cumprod(*args, **kwargs)

    @return_pyranges_if_possible
    def cumsum(self, *args, **kwargs):
        return super().cumsum(*args, **kwargs)

    @return_pyranges_if_possible
    def describe(self, *args, **kwargs):
        return super().describe(*args, **kwargs)

    @return_pyranges_if_possible
    def diff(self, *args, **kwargs):
        return super().diff(*args, **kwargs)

    @return_pyranges_if_possible
    def ewm(self, *args, **kwargs):
        return super().ewm(*args, **kwargs)

    @return_pyranges_if_possible
    def expanding(self, *args, **kwargs):
        return super().expanding(*args, **kwargs)

    @return_pyranges_if_possible
    def ffill(self, *args, **kwargs):
        return super().ffill(*args, **kwargs)

    @return_pyranges_if_possible
    def fillna(self, *args, **kwargs):
        return super().fillna(*args, **kwargs)

    @return_pyranges_if_possible
    def filter(self, *args, **kwargs):
        return super().filter(*args, **kwargs)

    @return_pyranges_if_possible
    def first(self, *args, **kwargs):
        return super().first(*args, **kwargs)

    @return_pyranges_if_possible
    def get_group(self, *args, **kwargs):
        return super().get_group(*args, **kwargs)

    @return_pyranges_if_possible
    def head(self, *args, **kwargs):
        return super().head(*args, **kwargs)

    @return_pyranges_if_possible
    def idxmax(self, *args, **kwargs):
        return super().idxmax(*args, **kwargs)

    @return_pyranges_if_possible
    def idxmin(self, *args, **kwargs):
        return super().idxmin(*args, **kwargs)

    @return_pyranges_if_possible
    def last(self, *args, **kwargs):
        return super().last(*args, **kwargs)

    @return_pyranges_if_possible
    def max(self, *args, **kwargs):
        return super().max(*args, **kwargs)

    @return_pyranges_if_possible
    def mean(self, *args, **kwargs):
        return super().mean(*args, **kwargs)

    @return_pyranges_if_possible
    def median(self, *args, **kwargs):
        return super().median(*args, **kwargs)

    @return_pyranges_if_possible
    def min(self, *args, **kwargs):
        return super().min(*args, **kwargs)

    @return_pyranges_if_possible
    def ngroup(self, *args, **kwargs):
        return super().ngroup(*args, **kwargs)

    @return_pyranges_if_possible
    def nunique(self, *args, **kwargs):
        return super().nunique(*args, **kwargs)

    @return_pyranges_if_possible
    def ohlc(self, *args, **kwargs):
        return super().ohlc(*args, **kwargs)

    @return_pyranges_if_possible
    def pct_change(self, *args, **kwargs):
        return super().pct_change(*args, **kwargs)

    @return_pyranges_if_possible
    def pipe(self, *args, **kwargs):
        return super().pipe(*args, **kwargs)

    @return_pyranges_if_possible
    def prod(self, *args, **kwargs):
        return super().prod(*args, **kwargs)

    @return_pyranges_if_possible
    def quantile(self, *args, **kwargs):
        return super().quantile(*args, **kwargs)

    @return_pyranges_if_possible
    def rank(self, *args, **kwargs):
        return super().rank(*args, **kwargs)

    @return_pyranges_if_possible
    def resample(self, *args, **kwargs):
        return super().resample(*args, **kwargs)

    @return_pyranges_if_possible
    def rolling(self, *args, **kwargs):
        return super().rolling(*args, **kwargs)

    @return_pyranges_if_possible
    def sample(self, *args, **kwargs):
        return super().sample(*args, **kwargs)

    @return_pyranges_if_possible
    def sem(self, *args, **kwargs):
        return super().sem(*args, **kwargs)

    @return_pyranges_if_possible
    def shift(self, *args, **kwargs):
        return super().shift(*args, **kwargs)

    @return_pyranges_if_possible
    def size(self, *args, **kwargs):
        return super().size(*args, **kwargs)

    @return_pyranges_if_possible
    def skew(self, *args, **kwargs):
        return super().skew(*args, **kwargs)

    @return_pyranges_if_possible
    def std(self, *args, **kwargs):
        return super().std(*args, **kwargs)

    @return_pyranges_if_possible
    def sum(self, *args, **kwargs):
        return super().sum(*args, **kwargs)

    @return_pyranges_if_possible
    def tail(self, *args, **kwargs):
        return super().tail(*args, **kwargs)

    @return_pyranges_if_possible
    def take(self, *args, **kwargs):
        return super().take(*args, **kwargs)

    @return_pyranges_if_possible
    def transform(self, *args, **kwargs):
        return super().transform(*args, **kwargs)

    @return_pyranges_if_possible
    def value_counts(self, *args, **kwargs):
        return super().value_counts(*args, **kwargs)

    @return_pyranges_if_possible
    def var(self, *args, **kwargs):
        return super().var(*args, **kwargs)

