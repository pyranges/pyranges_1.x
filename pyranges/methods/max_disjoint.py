import pandas as pd
from sorted_nearest import max_disjoint  # type: ignore[import]


def _max_disjoint(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    if df.empty:
        return df

    slack = kwargs.get("slack", 0)

    cdf = df.sort_values("End")
    # important: sorted_nearest interprets slack differently than pyranges
    # 0 slack in sorted_nearest means that bookended intervals are clustered
    # together, while in pyranges it means that they are not.
    idx = max_disjoint(cdf.index.values, cdf.Start.values, cdf.End.values, slack - 1)

    return cdf.reindex(idx)
