import pandas as pd
from sorted_nearest import annotate_clusters, cluster_by  # type: ignore[import]

from pyranges.names import END_COL, START_COL


def _cluster(
    df: pd.DataFrame,
    cluster_column: str,
    count_column: str | None = None,
    slack: int = 0,
    **_,
) -> pd.DataFrame:
    if df.empty:
        return None

    cdf = df.sort_values(START_COL)

    ids = annotate_clusters(cdf[START_COL].to_numpy(), cdf[END_COL].values, slack=slack)

    cdf.insert(df.shape[1], cluster_column, ids)

    if count_column:
        _count = cdf.groupby(cluster_column).Cluster.count()
        _count.name = count_column
        cdf = cdf.merge(_count, on=cluster_column)

    return cdf


def _cluster_by(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    if df.empty:
        return None

    slack = kwargs.get("slack", 0)
    count = kwargs.get("count_column", None)
    cluster_column = kwargs.get("cluster_column", False)

    by = kwargs["by"]

    cdf = df.sort_values([by, "Start"]) if isinstance(by, str) else df.sort_values([*by, "Start"])

    if isinstance(by, str):
        new_ids = (cdf[by] != cdf[by].shift()).cumsum()
    else:
        new_ids = (cdf[by] != cdf[by].shift()).any(axis=1).cumsum()

    cdf.insert(cdf.shape[1], "ClusterBy", new_ids)

    cdf = cdf.sort_values(["ClusterBy", "Start"])

    ids = cluster_by(cdf.Start.values, cdf.End.values, cdf.ClusterBy.values, slack)

    cdf = cdf.drop("ClusterBy", axis=1)
    cdf.insert(cdf.shape[1], cluster_column, ids)

    if count:
        _count = cdf.groupby("Cluster").Cluster.count()
        _count.name = "Count"
        cdf = cdf.merge_overlaps(_count, on="Cluster")

    return cdf
