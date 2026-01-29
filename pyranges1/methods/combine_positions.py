import numpy as np
import pandas as pd


def _intersect_interval_columns(
    starts: pd.Series,
    ends: pd.Series,
    starts2: pd.Series,
    ends2: pd.Series,
) -> tuple[pd.Series, pd.Series]:
    new_starts = pd.Series(
        np.where(starts > starts2.to_numpy(), starts, starts2),
        index=starts.index,
    )
    new_ends = pd.Series(
        np.where(ends < ends2.to_numpy(), ends, ends2),
        index=ends.index,
    )
    return new_starts, new_ends


def _union_interval_columns(
    starts: pd.Series,
    ends: pd.Series,
    starts2: pd.Series,
    ends2: pd.Series,
) -> tuple[pd.Series, pd.Series]:
    new_starts = pd.Series(
        np.where(starts < starts2.to_numpy(), starts, starts2),
        index=starts.index,
    )
    new_ends = pd.Series(
        np.where(ends > ends2.to_numpy(), ends, ends2),
        index=ends.index,
    )
    return new_starts, new_ends


def _swap_interval_columns(
    starts: pd.Series,  # noqa: ARG001
    ends: pd.Series,  # noqa: ARG001
    starts2: pd.Series,
    ends2: pd.Series,
) -> tuple[pd.Series, pd.Series]:
    return starts2, ends2
