import pandas as pd
from sorted_nearest import max_disjoint  # type: ignore[import]


def _max_disjoint(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    if df.empty:
        return df

    slack = kwargs.get("slack", 0)

    cdf = df.sort_values("End")

    idx = max_disjoint(cdf.index.values, cdf.Start.values, cdf.End.values, slack)

    return cdf.reindex(idx)
