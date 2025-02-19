import numpy as np
import pandas as pd
from natsort import natsorted


def sort_one_by_one(d: pd.DataFrame, col1: str, col2: str) -> pd.DataFrame:
    """Equivalent to pd.sort_values(by=[col1, col2]), but faster."""
    d = d.sort_values(by=[col2])
    return d.sort_values(by=[col1], kind="mergesort")


def sort_factorize_dict(df: pd.DataFrame, by: list[str], *, use_natsort=True) -> np.ndarray:
    """Return a numpy array of group-IDs for each row of df.

    The groups are defined by `by` columns and ordered by natsort.
    """
    # 1. Collect unique combinations of columns `by`.
    #    Each row in `unique_combos` is a unique combination.
    unique_combos = df[by].drop_duplicates()

    # 2. Convert each row to a tuple so natsort can compare them easily.
    combos_as_tuples = [tuple(x) for x in unique_combos.to_numpy()]

    # 3. Perform a natural sort on these tuples.
    combos_sorted = natsorted(combos_as_tuples) if use_natsort else sorted(combos_as_tuples)

    # 4. Build a dictionary mapping combination -> group_ID.
    combo_to_id = {combo: i for i, combo in enumerate(combos_sorted)}

    # 5. Map each row in `df` back to its group_ID.
    #    Convert each row to a tuple, then look up in the dictionary.
    #    One efficient way is to vectorize the tuple-conversion first,
    #    then map once.

    # Convert all df[by] rows to tuples
    all_as_tuples = df[by].to_numpy().tolist()  # list of lists
    all_as_tuples = list(map(tuple, all_as_tuples))

    # Map each tuple to an ID
    return np.array([combo_to_id[tup] for tup in all_as_tuples], dtype=np.uint32)
