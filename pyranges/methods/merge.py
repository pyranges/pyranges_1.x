import pandas as pd
from sorted_nearest import find_clusters  # type: ignore[import]

from pyranges.names import END_COL, START_COL


def _merge(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    if df.empty:
        return None

    slack = kwargs.get("slack", 0)
    by = kwargs["by"]

    cdf = df.sort_values(START_COL)

    starts, ends, number = find_clusters(cdf.Start.values, cdf.End.values, slack)

    by_values = df.head(1).squeeze()[by].to_dict()

    cluster_df = pd.DataFrame(
        {
            START_COL: starts,
            END_COL: ends,
        }
        | by_values,
    )
    # Sort columns in the original order of the dataframe.
    cluster_df = cluster_df[[c for c in cdf.columns if c in cluster_df]]

    if kwargs["count_col"]:
        cluster_df.insert(cluster_df.shape[1], kwargs["count_col"], number)

    return cluster_df
