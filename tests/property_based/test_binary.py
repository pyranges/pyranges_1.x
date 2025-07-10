import subprocess  # nosec
import tempfile
from io import StringIO
from itertools import product

import numpy as np
import pandas as pd
import pytest
from hypothesis import (
    HealthCheck,
    given,
    reproduce_failure,
    settings,
)

from tests.helpers import assert_df_equal
from tests.property_based.hypothesis_helper import (
    deadline,
    dfs_min,
    dfs_min2,
    max_examples,
)

strand_behavior = ["ignore", "same", "opposite"]
no_opposite = [False, "same"]


def run_bedtools(command, gr, gr2, strand_behavior, nearest_overlap=False, nearest_how=None, ties=""):
    bedtools_strand = {False: "", "same": "-s", "opposite": "-S"}[strand_behavior]
    bedtools_overlap = {True: "", False: "-io"}[nearest_overlap]
    bedtools_how = {"upstream": "-id", "downstream": "-iu", None: ""}[nearest_how] + " -D a"
    # print("bedtools how:", bedtools_how)
    ties = "-t " + ties if ties else ""

    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = f"{temp_dir}/f1.bed"
        f2 = f"{temp_dir}/f2.bed"
        gr.to_csv(f1, sep="\t", header=False, index=False)
        gr2.to_csv(f2, sep="\t", header=False, index=False)

        cmd = command.format(
            f1=f1,
            f2=f2,
            strand=bedtools_strand,
            overlap=bedtools_overlap,
            bedtools_how=bedtools_how,
            ties=ties,
        )
        # ignoring the below line in bandit as only strings created by
        # the test suite is run here; no user input ever sought
        return subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()  # nosec  # nosec


def read_bedtools_result_set_op(bedtools_result, strand_behavior):
    if strand_behavior:
        usecols = [0, 1, 2, 5]
        names = ["Chromosome", "Start", "End", "Strand"]
    else:
        usecols = [0, 1, 2]
        names = ["Chromosome", "Start", "End"]

    return pd.read_csv(
        StringIO(bedtools_result),
        header=None,
        usecols=usecols,
        names=names,
        # dtype={
        #     "Start": np.int32,
        #     "End": np.int32
        # },
        sep="\t",
    )


def compare_results(bedtools_df, result) -> None:
    # from pydbg import dbg
    # dbg(bedtools_df.dtypes)
    # dbg(result.df.dtypes)

    if not bedtools_df.empty:
        assert_df_equal(result, bedtools_df)
    else:
        assert bedtools_df.empty == result.empty


def compare_results_nearest(bedtools_df, result) -> None:
    if not bedtools_df.empty:
        bedtools_df = bedtools_df[bedtools_df.Distance != -1]

    result = result.df

    if len(result) != 0:
        bedtools_df = bedtools_df.sort_values(["Start", "End", "Distance"])
        result = result.sort_values(["Start", "End", "Distance"])
        result_df = result[["Chromosome", "Start", "End", "Strand", "Distance"]]
        assert_df_equal(result_df, bedtools_df)
    else:
        assert bedtools_df.empty


@pytest.mark.bedtools
@pytest.mark.parametrize("strand_behavior", no_opposite)
@settings(
    max_examples=max_examples,
    deadline=deadline,
    print_blob=True,
)
@given(gr=dfs_min(), gr2=dfs_min())  # pylint: disable=no-value-for-parameter
def test_set_intersect(gr, gr2, strand_behavior) -> None:
    set_intersect_command = "bedtools intersect {strand} -a <(sort -k1,1 -k2,2n {f1} | bedtools merge {strand} -c 4,5,6 -o first -i -) -b <(sort -k1,1 -k2,2n {f2} | bedtools merge {strand} -c 4,5,6 -o first -i -)"
    bedtools_result = run_bedtools(set_intersect_command, gr, gr2, strand_behavior)

    bedtools_df = read_bedtools_result_set_op(bedtools_result, strand_behavior)

    result = gr.set_intersect_overlaps(gr2, strand_behavior=strand_behavior)

    compare_results(bedtools_df, result)


@pytest.mark.bedtools
@pytest.mark.parametrize("strand_behavior", no_opposite)
@settings(
    max_examples=max_examples,
    deadline=deadline,
    print_blob=True,
)
@given(gr=dfs_min(), gr2=dfs_min())  # pylint: disable=no-value-for-parameter
# @reproduce_failure('4.15.0', b'AXicY2RAA4zoAgAAVQAD')
def test_set_union(gr, gr2, strand_behavior) -> None:
    set_union_command = "cat {f1} {f2} | bedtools sort | bedtools merge {strand} -c 4,5,6 -o first -i -"  # set_union_command = "bedtools merge {strand} -c 4,5,6 -o first -i {f1}"
    bedtools_result = run_bedtools(set_union_command, gr, gr2, strand_behavior)

    bedtools_df = read_bedtools_result_set_op(bedtools_result, strand_behavior)

    result = gr.set_union_overlaps(gr2, strand_behavior=strand_behavior)

    compare_results(bedtools_df, result)


@pytest.mark.bedtools
@pytest.mark.parametrize("strand_behavior", strand_behavior)
@settings(
    max_examples=max_examples,
    deadline=deadline,
    print_blob=True,
)
@given(gr=dfs_min(), gr2=dfs_min())  # pylint: disable=no-value-for-parameter
# @reproduce_failure('4.32.2', b'AXicY2RAA4wQzIgiCAAAgAAF')
# @reproduce_failure('5.5.4', b'AXicY2RABYyMEAqKGRgAAHMABg==')
@reproduce_failure("6.126.0", b"AXicY3ROKBJjYHBkhCEGRiCJzAcARwQDfg==")
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

    compare_results(bedtools_df, result)


@pytest.mark.bedtools
@pytest.mark.parametrize("strand_behavior", strand_behavior)
@settings(
    max_examples=max_examples,
    print_blob=True,
    deadline=deadline,
)
@given(gr=dfs_min(), gr2=dfs_min())  # pylint: disable=no-value-for-parameter
# @reproduce_failure('4.15.0', b'AXicY2RABoxghAoAAGkABA==')
# @reproduce_failure('4.15.0', b'AXicY2RABoxgxAAjQQAAAG8ABQ==')
# @reproduce_failure('4.15.0', b'AXicY2RABqwMDIwMaAAAALkACA==')
# @reproduce_failure('4.15.0', b'AXicY2RAA4xIJAgAAABcAAQ=')
# reproduce_failure('4.15.0', b'AXicY2RAAEYGhv9AkhHGgQIAFHQBBQ==')
# @reproduce_failure('4.15.0', b'AXicY2QAAUYGGGCEYIQAVAgAALUACA==')
def test_coverage(gr, gr2, strand_behavior) -> None:
    coverage_command = "bedtools coverage {strand} -a {f1} -b {f2}"

    bedtools_result = run_bedtools(coverage_command, gr, gr2, strand_behavior)

    bedtools_df = pd.read_csv(
        StringIO(bedtools_result),
        header=None,
        usecols=[0, 1, 2, 3, 4, 5, 6, 9],
        names=["Chromosome", "Start", "End", "Name", "Score", "Strand", "NumberOverlaps", "FractionOverlaps"],
        dtype={"FractionOverlap": np.float64},
        sep="\t",
    )

    result = gr.coverage(gr2, strand_behavior=strand_behavior)

    # assert len(result) > 0
    assert np.all(bedtools_df.NumberOverlaps.values == result.NumberOverlaps.values)
    np.testing.assert_allclose(bedtools_df.FractionOverlaps, result.FractionOverlaps, atol=1e-5)
    # compare_results(bedtools_df, result)


# @pytest.mark.bedtools
# @pytest.mark.parametrize("strand_behavior", strand_behavior)
# @settings(
#     max_examples=max_examples,
#     deadline=deadline,
#     suppress_health_check=HealthCheck.all())
# @given(gr=dfs_min(), gr2=dfs_min())  # pylint: disable=no-value-for-parameter
# @reproduce_failure('4.15.0', b'AXicY2RgYGAEIzgAsRkZUfkMDAAA2AAI')
# def test_no_intersect(gr, gr2, strand_behavior):

#     intersect_command = "bedtools intersect -v {strand} -a {f1} -b {f2}"

#     bedtools_result = run_bedtools(intersect_command, gr, gr2, strand_behavior)

#     bedtools_df = pd.read_csv(
#         StringIO(bedtools_result),
#         header=None,
#         names="Chromosome Start End Name Score Strand".split(),
#         sep="\t")

#     # bedtools bug: https://github.com/arq5x/bedtools2/issues/719
#     result = gr.no_overlap(gr2, strand_behavior=strand_behavior)

#     from pydbg import dbg
#     dbg(result)
#     dbg(bedtools_df)

#     # result2 = gr.intersect(gr2, strand_behavior)

#     compare_results(bedtools_df, result)


@pytest.mark.bedtools
@pytest.mark.parametrize("strand_behavior", ["same", "opposite", False])
@settings(
    max_examples=max_examples,
    deadline=deadline,
    print_blob=True,
)
@given(gr=dfs_min(), gr2=dfs_min())  # pylint: disable=no-value-for-parameter
# @reproduce_failure('4.5.7', b'AXicLYaJCQAACIS0/YfuuQRRAbVG94Dk5LHSBgJ3ABU=')
# @reproduce_failure('4.15.0', b'AXicY2QAAUYGGAVlIQAAAIIABQ==')
def test_subtraction(gr, gr2, strand_behavior) -> None:
    subtract_command = "bedtools subtract {strand} -a {f1} -b {f2}"

    bedtools_result = run_bedtools(subtract_command, gr, gr2, strand_behavior)

    bedtools_df = pd.read_csv(
        StringIO(bedtools_result),
        header=None,
        names=["Chromosome", "Start", "End", "Name", "Score", "Strand"],
        sep="\t",
    )

    result = gr.range_subtract(gr2, strand_behavior=strand_behavior)

    compare_results(bedtools_df, result)


nearest_hows = [None, "upstream", "downstream"]
overlaps = [True, False]


@pytest.mark.bedtools
@pytest.mark.parametrize(
    ("nearest_how", "overlap", "strand_behavior"), product(nearest_hows, overlaps, strand_behavior)
)
@settings(
    max_examples=max_examples,
    deadline=deadline,
    print_blob=True,
)
@given(gr=dfs_min(), gr2=dfs_min())  # pylint: disable=no-value-for-parameter
def test_nearest(gr, gr2, nearest_how, overlap, strand_behavior) -> None:
    nearest_command = "bedtools closest {bedtools_how} {strand} {overlap} -t first -d -a <(sort -k1,1 -k2,2n {f1}) -b <(sort -k1,1 -k2,2n {f2})"

    bedtools_result = run_bedtools(nearest_command, gr, gr2, strand_behavior, overlap, nearest_how)

    bedtools_df = pd.read_csv(
        StringIO(bedtools_result),
        header=None,
        names=["Chromosome", "Start", "End", "Strand", "Chromosome2", "Distance"],
        usecols=[0, 1, 2, 5, 6, 12],
        sep="\t",
    )

    bedtools_df.Distance = bedtools_df.Distance.abs()

    bedtools_df = bedtools_df[bedtools_df.Chromosome2 != "."]
    bedtools_df = bedtools_df.drop("Chromosome2", axis=1)

    result = gr.nearest_ranges(gr2, strand_behavior=strand_behavior, overlap=overlap, how=nearest_how)

    compare_results_nearest(bedtools_df, result)


@pytest.mark.bedtools
@pytest.mark.parametrize("strand_behavior", no_opposite)
@settings(
    max_examples=max_examples,
    deadline=deadline,
    print_blob=True,
    suppress_health_check=list(HealthCheck),
)
@given(gr=dfs_min(), gr2=dfs_min())  # pylint: disable=no-value-for-parameter
def test_jaccard(gr, gr2, strand_behavior) -> None:
    """Bedtools segfaults."""
    jaccard_command = (  # noqa: F841
        "bedtools jaccard {strand}  -a <(sort -k1,1 -k2,2n {f1}) -b <(sort -k1,1 -k2,2n {f2})"
    )

    #     # https://github.com/arq5x/bedtools2/issues/645
    #     # will make tests proper when bedtools is fixed
    result = gr.stats.jaccard(gr2, strand_behavior=strand_behavior)

    assert 0 <= result <= 1


@pytest.mark.bedtools
@pytest.mark.parametrize("strand_behavior", strand_behavior)
@settings(
    max_examples=max_examples,
    deadline=deadline,
    print_blob=True,
)
@given(gr=dfs_min(), gr2=dfs_min())  # pylint: disable=no-value-for-parameter
def test_join(gr, gr2, strand_behavior) -> None:
    join_command = "bedtools intersect {strand} -wo -a {f1} -b {f2}"

    bedtools_result = run_bedtools(join_command, gr, gr2, strand_behavior)

    bedtools_df = pd.read_csv(
        StringIO(bedtools_result),
        header=None,
        sep="\t",
        names=[
            "Chromosome",
            "Start",
            "End",
            "Name",
            "Score",
            "Strand",
            "Chromosome_b",
            "Start_b",
            "End_b",
            "Name_b",
            "Score_b",
            "Strand_b",
            "Overlap",
        ],
        dtype={"Chromosome": "category", "Strand": "category"},
    ).drop(["Chromosome_b", "Overlap"], axis=1)

    result = gr.join_overlaps(gr2, strand_behavior=strand_behavior)

    if result.df.empty:
        assert bedtools_df.empty
    else:
        assert_df_equal(result.df, bedtools_df)


@pytest.mark.bedtools
@settings(
    max_examples=max_examples,
    deadline=deadline,
    print_blob=True,
)
@given(gr=dfs_min2(), gr2=dfs_min2())  # pylint: disable=no-value-for-parameter
def test_reldist(gr, gr2) -> None:
    reldist_command = "bedtools reldist -a <(sort -k1,1 -k2,2n {f1}) -b <(sort -k1,1 -k2,2n {f2})"

    bedtools_result = run_bedtools(reldist_command, gr, gr2, False)
    bedtools_result = pd.read_csv(StringIO(bedtools_result), sep="\t")

    gr.stats.relative_distance(gr2)

    # bug in bedtools, therefore not testing this properly
    # https://github.com/arq5x/bedtools2/issues/711

    assert 1


new_pos = ["union"]  # ["intersection", "union"]


# @pytest.mark.parametrize("strand_behavior,new_pos", product(
#     strand_behavior, new_pos))
# @settings(
#     max_examples=max_examples,
#     deadline=deadline,
#     print_blob=True,
#     suppress_health_check=HealthCheck.all())
# @given(gr=dfs_min(), gr2=dfs_min())  # pylint: disable=no-value-for-parameter
# def test_join_new_pos(gr, gr2, strand_behavior, new_pos):

#     result = gr.join(gr2, strand_behavior=strand_behavior).new_position(new_pos)

#     import numpy as np
#     result2 = gr.join(gr2, strand_behavior=strand_behavior)

#     if result.df.empty:
#         assert result2.df.empty
#     else:
#         if new_pos == "union":
#             new_starts = np.minimum(result2.Start, result2.Start_b)
#             new_ends = np.maximum(result2.End, result2.End_b)
#         else:
#             new_starts = np.maximum(result2.Start, result2.Start_b)
#             new_ends = np.minimum(result2.End, result2.End_b)
#         assert list(result.Start.values) == list(new_starts)
#         assert list(result.End.values) == list(new_ends)


# @pytest.mark.parametrize("strand", [True, False])
# @settings(
#     max_examples=max_examples,
#     deadline=deadline,
#     suppress_health_check=HealthCheck.all())
# @given(gr=dfs_min_with_gene_id())  # pylint: disable=no-value-for-parameter
# def test_introns(gr, strand):

#     result = gr.features.introns()
#     print(result)

#     df = gr.df

#     grs = []
#     for g, gdf in df.groupby("ID"):
#         grs.append(pr.PyRanges(gdf))

#     expected = pr.concat([gr.merge() for gr in grs]).df

#     print(expected)
#     print(result)

#     assert_df_equal(result, expected)

k_nearest_ties = ["first", "last", None]
# k_nearest_ties = ["first", None]
k_nearest_ties = ["last"]

k_nearest_params = reversed(list(product(nearest_hows, [True, False], strand_behavior, k_nearest_ties)))


# @pytest.mark.bedtools
# @pytest.mark.explore
# @pytest.mark.parametrize("nearest_how,overlap,strand_behavior,ties", k_nearest_params)  #
# @settings(
#     max_examples=max_examples,
#     deadline=deadline,
#     print_blob=True,
#     suppress_health_check=HealthCheck.all(),
# )
# @given(gr=dfs_min(), gr2=dfs_min())  # pylint: disable=no-value-for-parameter
# # @reproduce_failure('4.43.5', b'AXicY2RAA4zoTAAAWwAE')
# def test_k_nearest(gr, gr2, nearest_how, overlap, strand_behavior, ties):
#     print("-----" * 20)
#
#     # gr = gr.apply(lambda df: df.astype({"Start": np.int32, "End": np.int32}))
#     # gr2 = gr2.apply(lambda df: df.astype({"Start": np.int32, "End": np.int32}))
#
#     # print(gr)
#     # print(gr2)
#
#     nearest_command = "bedtools closest -k 2 {bedtools_how} {strand} {overlap} {ties} -a <(sort -k1,1 -k2,2n {f1}) -b <(sort -k1,1 -k2,2n {f2})"
#
#     bedtools_result = run_bedtools(
#         nearest_command, gr, gr2, strand_behavior, overlap, nearest_how, ties
#     )
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
#     # cannot test with k > 1 because bedtools algo has different syntax
#     # cannot test keep_duplicates "all" or None/False properly, as the semantics is different for bedtools
#     result = gr.k_nearest(
#         gr2, k=2, strand_behavior=strand_behavior, overlap=overlap, how=nearest_how, ties=ties
#     )
#
#     # result = result.apply(lambda df: df.astype({"Start": np.int64, "End": np.int64, "Distance": np.int64}))
#     if len(result):
#         result.Distance = result.Distance.abs()
#     print("bedtools " * 5)
#     print(bedtools_df)
#     print("result " * 5)
#     print(result)
#
#     compare_results_nearest(bedtools_df, result)


# @settings(
#     max_examples=max_examples,
#     deadline=deadline,
#     print_blob=True,
#     suppress_health_check=HealthCheck.all())
# @given(gr=dfs_min())  # pylint: disable=no-value-for-parameter
# def test_k_nearest_nearest_self_same_size(gr):

#     result = gr.k_nearest(
#         gr, k=1, strand_behavior=None, overlap=True, how=None, ties="first")

#     assert len(result) == len(gr)


# @settings(
#     max_examples=max_examples,
#     deadline=deadline,
#     print_blob=True,
#     suppress_health_check=HealthCheck.all(),
# )
# @given(gr=dfs_min(), gr2=dfs_min())  # pylint: disable=no-value-for-parameter
# def test_k_nearest_1_vs_nearest(gr, gr2):
#     result_k = gr.k_nearest(gr2, k=1, strand_behavior=None, overlap=True, how=None)
#     if len(result_k) > 0:
#         result_k.Distance = result_k.Distance.abs()
#
#     result_n = gr.nearest_ranges(gr2, strand_behavior=None, overlap=True, how=None)
#
#     if len(result_k) == 0 and len(result_n) == 0:
#         pass
#     else:
#         assert (result_k.sort().Distance.abs() == result_n.sort().Distance).all()
