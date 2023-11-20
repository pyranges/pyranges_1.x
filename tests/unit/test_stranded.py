import numpy as np

import pyranges as pr
from pyranges.names import FORWARD_STRAND, REVERSE_STRAND

np.random.seed(0)


def test_stranded():
    cpg = pr.data.cpg()
    exons = pr.data.exons()

    j = cpg.interval_join(exons)

    assert j.strand_values_valid

    j.Strand = "."

    assert not j.strand_values_valid

    j.Strand = np.random.choice([REVERSE_STRAND, FORWARD_STRAND], size=len(j))

    assert j.strand_values_valid


def test_unstrand():
    exons = pr.data.exons()

    cpg = pr.data.cpg()
    print(exons)
    print(exons.columns)
    x = exons.remove_strand()
    print(x.strand_values_valid)
    print(x)
    print(x.columns)
    for _, df in x:
        assert not df.index.duplicated().sum()

    # x = pr.PyRanges(x.df.reset_index(drop=True))
    cpg.interval_join(x)
