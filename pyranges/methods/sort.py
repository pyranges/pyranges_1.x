import pandas as pd


def sort_one_by_one(d: pd.DataFrame, col1: str, col2: str) -> pd.DataFrame:
    """Equivalent to pd.sort_values(by=[col1, col2]), but faster."""
    d = d.sort_values(by=[col2])
    return d.sort_values(by=[col1], kind="mergesort")
