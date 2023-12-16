import pandas as pd

import pyranges as pr
from numpy import nan


def test_join_issue_4_right() -> None:
    import numpy as np

    chromsizes = pr.example_data.chromsizes
    query_regions = pr.gf.tile_genome(chromsizes, int(1e6))
    signal_data = pr.example_data.chipseq
    signal_data["Score"] = np.random.randint(0, 100, len(signal_data))

    query_regions.interval_join(signal_data)


def test_join_issue_8():
    gd = {'Chromosome': ['chr1', 'chr1', 'chr1', 'chr1'], 'Start': [157, 584, 731, 821], 'End': [257, 684, 831, 921], 'Strand': ['-', '-', '-', '-']}
    md = {'Chromosome': ['chr1', 'chr1', 'chr1', 'chr1'], 'Start': [316, 793, 889, 795], 'End': [416, 893, 989, 895], 'Strand': ['+', '+', '+', '-']}

    g = pr.PyRanges(gd)
    m = pr.PyRanges(md)

    j = m.interval_join(g)

    print(j)
    expected_result = pd.DataFrame({'Chromosome': [nan, nan, 'chr1', 'chr1'], 'Start': [nan, nan, 795.0, 795.0], 'End': [nan, nan, 895.0, 895.0], 'Strand': [nan, nan, '-', '-'], 'Chromosome_b': ['chr1', 'chr1', 'chr1', 'chr1'], 'Start_b': [157, 584, 731, 821], 'End_b': [257, 684, 831, 921], 'Strand_b': ['-', '-', '-', '-']})
    assert j.equals(expected_result)


def test_join_issue_8_right():
    gd = {'Chromosome': ['chr1', 'chr1', 'chr1', 'chr1'], 'Start': [157, 584, 731, 821], 'End': [257, 684, 831, 921], 'Strand': ['-', '-', '-', '-']}
    md = {'Chromosome': ['chr1', 'chr1', 'chr1', 'chr1'], 'Start': [316, 793, 889, 795], 'End': [416, 893, 989, 895], 'Strand': ['+', '+', '+', '-']}

    g = pr.PyRanges(gd)
    m = pr.PyRanges(md)

    j = m.interval_join(g, join_type="right")
    print(j)

    expected_result = pd.DataFrame({'Chromosome': [nan, nan, 'chr1', 'chr1'], 'Start': [nan, nan, 795.0, 795.0], 'End': [nan, nan, 895.0, 895.0], 'Strand': [nan, nan, '-', '-'], 'Chromosome_b': ['chr1', 'chr1', 'chr1', 'chr1'], 'Start_b': [157, 584, 731, 821], 'End_b': [257, 684, 831, 921], 'Strand_b': ['-', '-', '-', '-']})
    assert j.equals(expected_result)
