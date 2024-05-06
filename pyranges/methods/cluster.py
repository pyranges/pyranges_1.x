import pandas as pd
from sorted_nearest import annotate_clusters  # type: ignore[import]

from pyranges.core.names import END_COL, START_COL


def _cluster(
    df: pd.DataFrame,
    cluster_column: str,
    count_column: str | None = None,
    slack: int = 0,
    **_,
) -> pd.DataFrame:
    if df.empty:
        return df

    # important: sorted_nearest interprets slack differently than pyranges
    # 0 slack in sorted_nearest means that bookended intervals are clustered
    # together, while in pyranges it means that they are not.
    ids = annotate_clusters(df[START_COL].to_numpy(), df[END_COL].to_numpy(), slack=slack - 1)

    df.insert(df.shape[1], cluster_column, ids)

    if count_column:
        _count = df.groupby(cluster_column)["Cluster"].count()
        _count.name = count_column
        df = df.merge(_count, on=cluster_column)

    return df
