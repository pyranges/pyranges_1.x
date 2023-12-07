"""Test that overriden pandas methods work and return PyRanges objects."""

import pandas as pd
import pytest

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


@pytest.fixture
def gr():
    return pr.PyRanges(
        {
            "Chromosome": ["chr1"] * 2,
            "Start": [0, 10],
            "End": [40, 20],
            "Gene": ["DonkeyKong", "Mario"],
            "Val": [50, 30],
        }
    )


def test_loc_set():
    gr = pr.PyRanges({"Chromosome": ["chr1"], "Start": [0], "End": [40]})

    gr.loc[gr.Chromosome == "chr1", "Start"] = 1
    assert isinstance(gr, pr.PyRanges)
    assert gr.Start[0] == 1


def test_groupby_agg(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.agg("first"), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.agg("first"), pr.PyRanges)


def test_groupby_aggregate(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.aggregate("first"), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.aggregate("first"), pr.PyRanges)


def test_groupby_apply(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.apply(lambda df: df.head()), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.apply(lambda df: df.head()), pr.PyRanges)


def test_groupby_bfill(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.bfill(), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.bfill(), pr.PyRanges)


def test_groupby_ffill(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.ffill(), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.ffill(), pr.PyRanges)


def test_groupby_fillna(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.fillna(0), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.fillna(0), pr.PyRanges)


def test_groupby_filter(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.filter(lambda df: True), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.filter(lambda df: True), pr.PyRanges)


def test_groupby_first(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.first(), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.first(), pr.PyRanges)


def test_groupby_get_group(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.get_group("chr1"), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.get_group("Mario"), pr.PyRanges)


def test_groupby_head(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.head(), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.head(), pr.PyRanges)


def test_groupby_last(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.last(), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.last(), pr.PyRanges)


def test_groupby_min(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.min(), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.min(), pr.PyRanges)


def test_groupby_pipe(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.pipe(lambda df: df).agg("first"), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.pipe(lambda df: df).agg("first"), pr.PyRanges)


def test_groupby_prod(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.prod("first"), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.prod(), pr.PyRanges)


def test_groupby_sample(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.sample(), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.sample(), pr.PyRanges)


def test_groupby_sum(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.sum(), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.sum(), pr.PyRanges)


def test_groupby_tail(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.tail(1), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.tail(1), pr.PyRanges)


def test_groupby_transform(gr):
    g = gr.groupby("Chromosome")
    assert isinstance(g.transform("first"), pd.DataFrame)

    g = gr.groupby("Gene")
    assert isinstance(g.transform("first"), pr.PyRanges)
