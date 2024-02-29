import pandas as pd

import pyranges as pr
import contextlib
import os


def test_tostring():
    with contextlib.redirect_stdout(open(os.devnull, "w+")):
        print(pr.example_data.f1)


def test_tostring_index_and_column_share_name():
    gr = pr.example_data.f1
    gr = gr.set_index("Chromosome", drop=False)
    with contextlib.redirect_stdout(open(os.devnull, "w+")):
        print(gr)
