import contextlib
import os

import pyranges as pr


def test_tostring() -> None:
    with contextlib.redirect_stdout(open(os.devnull, "w+")):
        pass


def test_tostring_index_and_column_share_name() -> None:
    gr = pr.example_data.f1
    gr = gr.set_index("Chromosome", drop=False)
    with contextlib.redirect_stdout(open(os.devnull, "w+")):
        pass
