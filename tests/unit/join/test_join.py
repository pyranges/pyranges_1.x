import numpy as np
from pandas.testing import assert_series_equal

import pyranges as pr


# with slack
def test_join_with_slack():
    gr1 = pr.from_args(chromosomes="chr1", starts=[0], ends=[10], strands="+")
    gr2 = pr.from_args(chromosomes="chr1", starts=[15], ends=[20], strands="+")

    result = gr1.interval_join(gr2, slack=10)
    assert not result.empty


def test_join_without_reordering():
    f1 = pr.PyRanges(
        {
            "Chromosome": ["chr1", "chr1", "chr1"],
            "Start": [3, 8, 5],
            "End": [6, 9, 7],
            "Name": ["interval1", "interval3", "interval2"],
        }
    )
    f2 = pr.PyRanges(
        {
            "Chromosome": ["chr1", "chr1"],
            "Start": [1, 6],
            "End": [2, 7],
            "Name": ["a", "b"],
        }
    )

    lj = f1.interval_join(f2, join_type="left")
    assert_series_equal(lj.Name, f1.Name)

    rj = f1.interval_join(f2, join_type="right")
    assert_series_equal(rj.Name_b, f2.Name, check_names=False)

    oj = f1.interval_join(f2, join_type="outer")
    assert list(oj.Name.fillna(-1)) == ["interval2", "interval1", "interval3", -1]
    assert list(oj.Name_b.fillna(-1)) == ["b", -1, -1, "a"]
