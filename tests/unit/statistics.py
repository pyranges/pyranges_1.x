import pyranges as pr


def test_mcc():
    gr1, gr2 = pr.example_data.chipseq(), pr.example_data.chipseq_background()
    g = pr.example_data.chromsizes()
    mcc = pr.stats.mcc([gr1, gr2], genome=g, labels=["chip", "bg"], strand=True)

    print(mcc)
