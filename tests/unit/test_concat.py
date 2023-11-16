#!/usr/bin/env python3

import pyranges as pr


def assert_equal_length_before_after(gr1, gr2):
    print("in test")
    l1 = len(gr1)
    l2 = len(gr2)
    c = pr.concat([gr1, gr2])

    if not gr1.stranded or not gr2.stranded:
        assert not c.stranded

    lc = len(c)
    assert l1 + l2 == lc


def test_concat_stranded_tranded(f1, f2):
    assert_equal_length_before_after(f1, f2)


def test_concat_unstranded_unstranded(f1, f2):
    assert_equal_length_before_after(f1.remove_strand(), f2.remove_strand())


def test_concat_stranded_unstranded(f1, f2):
    assert_equal_length_before_after(f1, f2.remove_strand())


def test_concat_unstranded_stranded(f1, f2):
    assert_equal_length_before_after(f1.remove_strand(), f2)
