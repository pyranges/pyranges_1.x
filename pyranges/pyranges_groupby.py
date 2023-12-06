import pandas as pd

from pyranges import GENOME_LOC_COLS


def return_pyranges_if_possible(method):
    def wrapper(*args, **kwargs):
        # Call the original groupby method
        result = method(*args, **kwargs)

        # Check if the result should be a MySpecialDataFrame
        if isinstance(result, pd.DataFrame) and set(GENOME_LOC_COLS).issubset(result.columns):
            return PyRangesGroupBy(result)
        else:
            return result

    return wrapper


class PyRangesGroupBy(pd.core.groupby.DataFrameGroupBy):

    @return_pyranges_if_possible
    def aggregate(self, func, *args, **kwargs):
        return super().aggregate(func, *args, **kwargs)

    @return_pyranges_if_possible
    def apply(self, func, *args, **kwargs):
        return super().apply(func, *args, **kwargs)

    # Apply the decorator to other methods as needed
