import pyranges as pr


def test_unstrand():
    gr = pr.data.chipseq()
    u = gr.remove_strand()
    assert not u.strand_values_valid
