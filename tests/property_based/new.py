import subprocess
import tempfile
from io import StringIO

import hypothesis.strategies as st
import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis.extra.pandas import column, data_frames, indexes

from pyranges import PyRanges
from pyranges.core.names import VALID_GENOMIC_STRAND_INFO

max_examples = 15
slow_max_examples = 5
deadline = None

starts = st.integers(min_value=1, max_value=int(1e7))
lengths = st.integers(min_value=1, max_value=int(1e4))

strands = st.sampled_from(VALID_GENOMIC_STRAND_INFO)

chromosomes = st.sampled_from([f"chr{e!s}" for e in [*list(range(1, 23)), "X", "Y", "M"]])


dfs_zero_length_not_allowed = data_frames(
    index=indexes(dtype=np.int64, min_size=1, unique=True, elements=starts),
    columns=[
        column("Chromosome", chromosomes),
        column("Start", elements=starts),
        column("End", elements=lengths),
        column("Strand", strands),
    ],
)


@st.composite
def nonempty_pyranges(draw):  # nosec
    df = draw(dfs_zero_length_not_allowed)
    df.loc[:, "End"] += df.Start

    df.insert(3, "Name", "a")
    df.insert(4, "Score", 0)

    return PyRanges(df)


def run_bedtools(command, gr, gr2, strand_behavior, nearest_overlap=False, nearest_how=None, ties=""):
    bedtools_strand_behavior = {"ignore": "", "same": "-s", "opposite": "-S"}[strand_behavior]
    bedtools_overlap = {True: "", False: "-io"}[nearest_overlap]
    bedtools_how = {"upstream": "-id", "downstream": "-iu", None: ""}[nearest_how] + " -D a"
    ties = "-t " + ties if ties else ""

    # file1_proc_sub = f"""<(echo -e '{gr.to_csv(sep="\t", header=False, index=False)}')"""
    # file2_proc_sub = f"""<(echo -e '{gr2.to_csv(sep="\t", header=False, index=False)}')"""

    # cmd = command.format(
    #     f1=file1_proc_sub,
    #     f2=file2_proc_sub,
    #     strand=bedtools_strand_behavior,
    #     overlap=bedtools_overlap,
    #     bedtools_how=bedtools_how,
    #     ties=ties,
    # )
    # return subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()  # nosec  # nosec

    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = f"{temp_dir}/f1.bed"
        f2 = f"{temp_dir}/f2.bed"
        gr.to_csv(f1, sep="\t", header=False, index=False)
        gr2.to_csv(f2, sep="\t", header=False, index=False)

        cmd = command.format(
            f1=f1,
            f2=f2,
            strand=bedtools_strand_behavior,
            overlap=bedtools_overlap,
            bedtools_how=bedtools_how,
            ties=ties,
        )
        # ignoring the below line in bandit as only strings created by
        # the test suite is run here; no user input ever sought
        return subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()  # nosec  # nosec


def assert_equal(result, bedtools_df) -> None:
    if result.empty and bedtools_df.empty:
        return
    result = PyRanges(result).sort_ranges().reset_index(drop=True)
    bedtools_df = PyRanges(bedtools_df).sort_ranges().reset_index(drop=True)
    pd.testing.assert_frame_equal(result, bedtools_df, check_exact=False, atol=1e-5)


strand_behavior = ["ignore", "same", "opposite"]


@pytest.mark.bedtools
@pytest.mark.parametrize("strand_behavior", strand_behavior)
@settings(
    max_examples=max_examples,
    print_blob=True,
)
@given(gr=nonempty_pyranges(), gr2=nonempty_pyranges())  # pylint: disable=no-value-for-parameter
def test_overlap(gr, gr2, strand_behavior) -> None:
    overlap_command = "bedtools intersect -u {strand} -a {f1} -b {f2}"

    bedtools_result = run_bedtools(overlap_command, gr, gr2, strand_behavior)

    bedtools_df = pd.read_csv(
        StringIO(bedtools_result),
        header=None,
        names=["Chromosome", "Start", "End", "Name", "Score", "Strand"],
        sep="\t",
    )

    result = gr.overlap(gr2, strand_behavior=strand_behavior)

    result = result.reset_index(drop=True)
    if result.empty and bedtools_df.empty:
        return
    assert_equal(result, bedtools_df)


@pytest.mark.bedtools
@pytest.mark.parametrize("strand_behavior", strand_behavior)
@settings(
    max_examples=max_examples,
    print_blob=True,
)
@given(gr=nonempty_pyranges(), gr2=nonempty_pyranges())  # pylint: disable=no-value-for-parameter
def test_coverage(gr, gr2, strand_behavior) -> None:
    coverage_command = "bedtools coverage {strand} -a {f1} -b {f2}"

    bedtools_result = run_bedtools(coverage_command, gr, gr2, strand_behavior)

    bedtools_df = PyRanges(
        pd.read_csv(
            StringIO(bedtools_result),
            header=None,
            usecols=[0, 1, 2, 3, 4, 5, 6, 9],
            names=["Chromosome", "Start", "End", "Name", "Score", "Strand", "NumberOverlaps", "FractionOverlaps"],
            dtype={"FractionOverlap": np.float64},
            sep="\t",
        )
    )

    result = gr.coverage(gr2, strand_behavior=strand_behavior)

    assert_equal(result, bedtools_df)


STRAND_BEHAVIOR_NO_OPPOSITE = ["ignore", "same"]


@pytest.mark.bedtools
@pytest.mark.parametrize("strand_behavior", STRAND_BEHAVIOR_NO_OPPOSITE)
@settings(
    max_examples=max_examples,
    print_blob=True,
)
@given(gr=nonempty_pyranges(), gr2=nonempty_pyranges())  # pylint: disable=no-value-for-parameter
def test_set_intersect(gr, gr2, strand_behavior) -> None:
    set_intersect_command = "bedtools intersect {strand} -a <(sort -k1,1 -k2,2n {f1} | bedtools merge {strand} -c 4,5,6 -o first -i -) -b <(sort -k1,1 -k2,2n {f2} | bedtools merge {strand} -c 4,5,6 -o first -i -)"
    bedtools_result = run_bedtools(set_intersect_command, gr, gr2, strand_behavior)

    bedtools_df = read_bedtools_result_set_op(bedtools_result, strand_behavior=strand_behavior)

    result = gr.set_intersect_overlaps(gr2, strand_behavior=strand_behavior)

    if result.empty and bedtools_df.empty:
        return

    assert_equal(result, bedtools_df)


@pytest.mark.bedtools
@pytest.mark.parametrize("strand_behavior", STRAND_BEHAVIOR_NO_OPPOSITE)
@settings(
    max_examples=max_examples,
    deadline=deadline,
    print_blob=True,
)
@given(gr=nonempty_pyranges(), gr2=nonempty_pyranges())
def test_set_union(gr, gr2, strand_behavior) -> None:
    set_union_command = "cat {f1} {f2} | bedtools sort | bedtools merge {strand} -c 4,5,6 -o first -i -"  # set_union_command = "bedtools merge {strand} -c 4,5,6 -o first -i {f1}"
    bedtools_result = run_bedtools(set_union_command, gr, gr2, strand_behavior)

    bedtools_df = read_bedtools_result_set_op(bedtools_result, strand_behavior)

    result = gr.set_union_overlaps(gr2, strand_behavior=strand_behavior)

    assert_equal(bedtools_df, result)


# @pytest.mark.bedtools
# @pytest.mark.parametrize(("nearest_how", "overlap", "strand_behavior"), product(nearest_hows, overlaps, strand_behavior))
# @settings(
#     max_examples=max_examples,
#     deadline=deadline,
#     print_blob=True,
# )
# @given(gr=nonempty_pyranges(), gr2=nonempty_pyranges())  # pylint: disable=no-value-for-parameter
# def test_nearest(gr, gr2, nearest_how, overlap, strand_behavior) -> None:
#     nearest_command = "bedtools closest {bedtools_how} {strand} {overlap} -t first -d -a <(sort -k1,1 -k2,2n {f1}) -b <(sort -k1,1 -k2,2n {f2})"
#
#     bedtools_result = run_bedtools(nearest_command, gr, gr2, strand_behavior, overlap, nearest_how)
#
#     bedtools_df = pd.read_csv(
#         StringIO(bedtools_result),
#         header=None,
#         names="Chromosome Start End Strand Chromosome2 Distance".split(),
#         usecols=[0, 1, 2, 5, 6, 12],
#         sep="\t",
#     )
#
#     bedtools_df.Distance = bedtools_df.Distance.abs()
#
#     bedtools_df = bedtools_df[bedtools_df.Chromosome2 != "."]
#     bedtools_df = bedtools_df.drop("Chromosome2", axis=1)
#
#     result = gr.nearest_ranges(gr2, strand_behavior=strand_behavior, overlap=overlap, how=nearest_how)
#
#     print("bedtools " * 5)
#     print(bedtools_df)
#     print("result " * 5)
#     print(result)
#
#     compare_results_nearest(bedtools_df, result)


@pytest.mark.bedtools
@pytest.mark.parametrize("strand_behavior", strand_behavior)
@settings(
    max_examples=max_examples,
    print_blob=True,
)
@given(gr=nonempty_pyranges(), gr2=nonempty_pyranges())  # pylint: disable=no-value-for-parameter
def test_subtraction(gr, gr2, strand_behavior) -> None:
    subtract_command = "bedtools subtract {strand} -a {f1} -b {f2}"

    bedtools_result = run_bedtools(subtract_command, gr, gr2, strand_behavior)

    bedtools_df = pd.read_csv(
        StringIO(bedtools_result),
        header=None,
        names=["Chromosome", "Start", "End", "Name", "Score", "Strand"],
        sep="\t",
    )

    result = gr.subtract_overlaps(gr2, strand_behavior=strand_behavior)

    assert_equal(result, bedtools_df)


def read_bedtools_result_set_op(bedtools_result, strand_behavior):
    if strand_behavior == "same":
        usecols = [0, 1, 2, 5]
        names = ["Chromosome", "Start", "End", "Strand"]
    elif strand_behavior == "ignore":
        usecols = [0, 1, 2]
        names = ["Chromosome", "Start", "End"]
    else:
        raise ValueError("Invalid strand behavior: " + strand_behavior)

    return pd.read_csv(
        StringIO(bedtools_result),
        header=None,
        usecols=usecols,
        names=names,
        sep="\t",
    )
