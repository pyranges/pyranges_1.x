import numpy as np

import pyranges as pr


def test_i32() -> None:
    gr = pr.example_data.aorta

    gr = gr.astype({"Start": np.int32, "End": np.int32})

    gr.merge_overlaps()
    assert 0
