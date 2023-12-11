import pandas as pd
from sorted_nearest import annotate_clusters  # type: ignore[import]

from pyranges.names import END_COL, START_COL


def _cluster(
    df: pd.DataFrame,
    cluster_column: str,
    count_column: str | None = None,
    slack: int = 0,
    **_,
) -> pd.DataFrame:
    if df.empty:
        return df

    cdf = df.sort_values(START_COL)

    ids = annotate_clusters(cdf[START_COL].to_numpy(), cdf[END_COL].values, slack=slack)

    cdf.insert(df.shape[1], cluster_column, ids)

    if count_column:
        _count = cdf.groupby(cluster_column).Cluster.count()
        _count.name = count_column
        cdf = cdf.merge(_count, on=cluster_column)

    return cdf
