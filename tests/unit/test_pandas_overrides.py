"""Test that overriden pandas methods work and return PyRanges objects."""

import pandas as pd

import pyranges as pr


def test_getitem():
    gr = pr.PyRanges({"Chromosome": ["chr1"], "Start": [0], "End": [40]})
    res = gr[gr.Chromosome == "chr1"]
    assert isinstance(res, pr.PyRanges)


def test_loc_get():
    gr = pr.PyRanges({"Chromosome": ["chr1"], "Start": [0], "End": [40]})
    res = gr.loc[gr.Chromosome == "chr1", ["Chromosome", "Start", "End"]]
    assert isinstance(res, pr.PyRanges)

    res = gr.loc[gr.Chromosome == "chr1", ["Chromosome", "Start"]]
    assert isinstance(res, pd.DataFrame)

    res = gr.loc[gr.Chromosome != "chr1"]
    assert isinstance(res, pr.PyRanges)


def test_loc_set():
    gr = pr.PyRanges({"Chromosome": ["chr1"], "Start": [0], "End": [40]})

    gr.loc[gr.Chromosome == "chr1", "Start"] = 1
    assert isinstance(gr, pr.PyRanges)
    assert gr.Start[0] == 1


def test_groupby():
    gr = pr.PyRanges({"Chromosome": ["chr1", "chr2"], "Start": [0, 4], "End": [40, 10]})
    res = gr.groupby("Chromosome").agg({"Chromosome": "first", "Start": "min", "End": "max"})
    print(res)

    assert 0



