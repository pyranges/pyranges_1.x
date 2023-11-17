import numpy as np

import pyranges as pr

np.random.seed(0)


def test_stranded():
    cpg = pr.data.cpg()
    exons = pr.data.exons()

    j = cpg.join_overlaps(exons)

    assert j.valid_strand

    j.Strand = "."

    assert not j.valid_strand

    j.Strand = np.random.choice("+ -".split(), size=len(j))

    assert j.valid_strand

    for _, df in j:
        assert len(df.Strand.drop_duplicates()) == 1


def test_unstrand():
    exons = pr.data.exons()

    cpg = pr.data.cpg()
    print(exons)
    print(exons.columns)
    x = exons.remove_strand()
    print(x.valid_strand)
    print(x)
    print(x.columns)
    for _, df in x:
        assert not df.index.duplicated().sum()

    # x = pr.PyRanges(x.df.reset_index(drop=True))
    cpg.join_overlaps(x)
