import numpy as np
import pandas as pd
import pytest

import pyranges1 as pr


def _base_frame_with_dtype(dtype: np.dtype) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Chromosome": pd.Series(["chr1", "chr1", "chr1", "chr1"], dtype="category"),
            "Start": pd.Series([1, 5, 20, 21], dtype=dtype),
            "End": pd.Series([10, 8, 25, 30], dtype=dtype),
            "Strand": pd.Series(["+", "+", "+", "+"], dtype="category"),
        },
    )


def test_cluster_overlaps_supports_uint32_issue_case() -> None:
    # Regression test for the historical unsupported dtype pair:
    # group dtype uint32 + position dtype uint32.
    gr = pr.PyRanges(_base_frame_with_dtype(np.uint32))

    out = gr.cluster_overlaps()

    assert out["Start"].dtype == np.dtype(np.uint32)
    assert out["End"].dtype == np.dtype(np.uint32)
    assert out["Cluster"].dtype == np.dtype(np.uint32)
    assert out["Cluster"].tolist() == [0, 0, 1, 1]


@pytest.mark.parametrize("dtype", [np.int8, np.uint8, np.int16, np.uint16, np.int32, np.uint32, np.int64])
def test_core_ops_support_common_integer_position_dtypes(dtype: np.dtype) -> None:
    gr = pr.PyRanges(_base_frame_with_dtype(dtype))

    clustered = gr.cluster_overlaps()
    merged = gr.merge_overlaps()

    assert clustered["Start"].dtype == np.dtype(dtype)
    assert clustered["End"].dtype == np.dtype(dtype)
    assert merged["Start"].dtype == np.dtype(dtype)
    assert merged["End"].dtype == np.dtype(dtype)
    assert clustered["Cluster"].dtype == np.dtype(np.uint32)


def test_uint64_outside_int64_range_raises() -> None:
    starts = np.array([2**63 + 1, 2**63 + 10], dtype=np.uint64)
    ends = starts + np.uint64(1)

    gr = pr.PyRanges(
        pd.DataFrame(
            {
                "Chromosome": ["chr1", "chr1"],
                "Start": starts,
                "End": ends,
            },
        ),
    )

    with pytest.raises(OverflowError, match="int64 range"):
        gr.cluster_overlaps()
