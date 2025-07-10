import subprocess  # nosec
import tempfile
from io import StringIO

import numpy as np
import pandas as pd
import pytest
from hypothesis import (
    HealthCheck,
    given,
    reproduce_failure,  # noqa: F401
    settings,
)
from natsort import natsorted  # type: ignore

import pyranges as pr
from tests.helpers import assert_df_equal
from tests.property_based.hypothesis_helper import (
    deadline,
    df_data,
    dfs_min,
    dfs_min_with_id,
    max_examples,
    selector,
)

# if environ.get("TRAVIS"):
#     max_examples = 100
#     deadline = None
# else:
#     max_examples = 1000
#     deadline = None

merge_command = "bedtools merge -o first,count -c 6,1 {} -i <(sort -k1,1 -k2,2n {})"


@pytest.mark.bedtools
@pytest.mark.parametrize("strand", [True, False])
@settings(
    max_examples=max_examples,
    deadline=deadline,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=dfs_min())  # pylint: disable=no-value-for-parameter
def test_merge(gr, strand) -> None:
    bedtools_strand = {True: "-s", False: ""}[strand]

    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = f"{temp_dir}/f1.bed"
        gr.df.to_csv(f1, sep="\t", header=False, index=False)

        cmd = merge_command.format(bedtools_strand, f1)

        # ignoring bandit security warning. All strings created by test suite
        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()  # nosec  # nosec

        if not strand:
            bedtools_df = pd.read_csv(
                StringIO(result),
                sep="\t",
                header=None,
                usecols=[0, 1, 2, 4],
                names=["Chromosome", "Start", "End", "Count"],
                dtype={"Chromosome": "category"},
            )
        else:
            bedtools_df = pd.read_csv(
                StringIO(result),
                sep="\t",
                header=None,
                names=["Chromosome", "Start", "End", "Strand", "Count"],
                dtype={"Chromosome": "category"},
            )

    result = gr.merge_overlaps(use_strand=strand)

    if not bedtools_df.empty:
        # need to sort because bedtools sometimes gives the result in non-natsorted chromosome order!
        if result.strand_valid:
            assert_df_equal(
                result.df.sort_values(["Chromosome", "Start", "Strand"]),
                bedtools_df.sort_values(["Chromosome", "Start", "Strand"]),
            )
        else:
            assert_df_equal(
                result.df.sort_values(["Chromosome", "Start"]),
                bedtools_df.sort_values(["Chromosome", "Start"]),
            )
    else:
        assert bedtools_df.empty == result.df.empty


cluster_command = "bedtools cluster {} -i <(sort -k1,1 -k2,2n {})"


@pytest.mark.bedtools
@pytest.mark.parametrize("strand", [True, False])
@settings(
    max_examples=max_examples,
    deadline=deadline,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=dfs_min())  # pylint: disable=no-value-for-parameter
def test_cluster(gr, strand) -> None:
    bedtools_strand = {True: "-s", False: ""}[strand]

    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = f"{temp_dir}/f1.bed"
        gr.df.to_csv(f1, sep="\t", header=False, index=False)

        cmd = cluster_command.format(bedtools_strand, f1)

        # ignoring bandit security warning. All strings created by test suite
        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()  # nosec  # nosec

        bedtools_df = pd.read_csv(
            StringIO(result),
            sep="\t",
            header=None,
            names=["Chromosome", "Start", "End", "Name", "Score", "Strand", "Cluster"],
            dtype={"Chromosome": "category"},
        )

    result = gr.cluster_overlaps(use_strand=strand)

    if not bedtools_df.empty:
        # need to sort because bedtools sometimes gives the result in non-natsorted chromosome order!
        sort_values = ["Chromosome", "Start", "Strand"] if result.strand_valid else ["Chromosome", "Start"]

        result_df = result.df.sort_values(sort_values)
        bedtools_df = bedtools_df.sort_values(sort_values)

        cluster_ids = dict(
            zip(
                result_df.Cluster.drop_duplicates(),
                bedtools_df.Cluster.drop_duplicates(),
                strict=False,
            ),
        )

        # bedtools gives different cluster ids than pyranges
        result_df.Cluster = result_df.Cluster.replace(cluster_ids)

        bedtools_df.Cluster = bedtools_df.Cluster.astype("int32")
        assert_df_equal(result_df.drop("Cluster", axis=1), bedtools_df.drop("Cluster", axis=1))
    else:
        assert bedtools_df.empty == result.df.empty


@pytest.mark.parametrize("strand", [True, False])
@settings(
    max_examples=max_examples,
    deadline=deadline,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=dfs_min_with_id())  # pylint: disable=no-value-for-parameter
def test_cluster_by(gr, strand) -> None:
    result = gr.cluster_overlaps(use_strand=strand, by="ID").df
    df = gr.df

    groupby = ["Chromosome", "Strand", "ID"] if strand else ["Chromosome", "ID"]

    grs = []

    for _, gdf in natsorted(df.groupby(groupby)):
        grs.append(pr.PyRanges(gdf))

    clusters = [gr.cluster_overlaps(use_strand=strand) for gr in grs]
    i = 1
    new_clusters = []
    for c in clusters:
        c.Cluster = i
        i += 1
        new_clusters.append(c)

    expected = pr.concat(new_clusters).df
    expected.loc[:, "Cluster"] = expected.Cluster.astype(np.int32)
    # expected = expected.drop_duplicates()

    assert_df_equal(result.drop("Cluster", axis=1), expected.drop("Cluster", axis=1))


@pytest.mark.parametrize("strand", [True, False])
@settings(
    max_examples=max_examples,
    deadline=deadline,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=dfs_min_with_id())  # pylint: disable=no-value-for-parameter
def test_merge_by(gr, strand) -> None:
    result = gr.merge_overlaps(by="ID").df.drop("ID", axis=1)

    df = gr.df

    grs = []
    for _, gdf in df.groupby("ID"):
        grs.append(pr.PyRanges(gdf))

    expected = pr.concat([gr.merge_overlaps() for gr in grs]).df

    assert_df_equal(result, expected)


makewindows_command = "bedtools makewindows -w 10 -b <(sort -k1,1 -k2,2n {})"


@pytest.mark.bedtools
@settings(
    max_examples=max_examples,
    print_blob=True,
    deadline=deadline,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=dfs_min())  # pylint: disable=no-value-for-parameter
# @reproduce_failure('5.5.4', b'AXicY2RgYGAEIzgAsRkBAFsABg==')
def test_windows(gr) -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        f1 = f"{temp_dir}/f1.bed"
        gr.df.to_csv(f1, sep="\t", header=False, index=False)

        cmd = makewindows_command.format(f1)

        # ignoring bandit security warning. All strings created by test suite
        result = subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()  # nosec  # nosec

        bedtools_df = pd.read_csv(
            StringIO(result),
            sep="\t",
            header=None,
            names=["Chromosome", "Start", "End"],
            dtype={"Chromosome": "category"},
        )

    result = gr.window_ranges(10)[["Chromosome", "Start", "End"]].remove_strand()

    if not bedtools_df.empty:
        assert_df_equal(result.df, bedtools_df)
    else:
        assert bedtools_df.empty == result.df.empty


@pytest.mark.parametrize("strand", [True, False])
@settings(
    max_examples=max_examples,
    deadline=deadline,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=df_data())  # pylint: disable=no-value-for-parameter
def test_init(gr, strand) -> None:
    c, s, e, strands = gr

    if strand:
        pr.PyRanges(chromosomes=c, starts=s, ends=e, strands=strands)
    else:
        pr.PyRanges(chromosomes=c, starts=s, ends=e)


chipseq = pr.example_data.chipseq()


@settings(
    max_examples=max_examples,
    deadline=deadline,
    suppress_health_check=HealthCheck.all(),
)
@given(selector=selector())  # pylint: disable=no-value-for-parameter
def test_getitem(selector) -> None:
    if len(selector) == 3:
        a, b, c = selector
        chipseq[a, b, c]
    elif len(selector) == 2:
        a, b = selector
        chipseq[a, b]
    elif len(selector) == 1:
        a = selector[0]
        chipseq[a]
    elif len(selector) == 0:
        pass
    else:
        msg = "Should never happen"
        raise Exception(msg)


@pytest.mark.bedtools
@settings(
    max_examples=max_examples,
    deadline=deadline,
    print_blob=True,
    suppress_health_check=HealthCheck.all(),
)
@given(gr=dfs_min())  # pylint: disable=no-value-for-parameter
def test_summary(gr) -> None:
    # merely testing that it does not error
    # contents are just (pandas) dataframe.describe()
    gr.summary()
