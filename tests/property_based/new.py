import subprocess
import tempfile

import hypothesis.strategies as st
import numpy as np
import pandas as pd
import pytest
from hypothesis import given, reproduce_failure, settings
from hypothesis.extra.pandas import column, data_frames, indexes

from pyranges import PyRanges
from pyranges.names import VALID_GENOMIC_STRAND_INFO
from tests.property_based.test_binary import read_bedtools_result_set_op

max_examples = 15
slow_max_examples = 5
deadline = None

starts = st.integers(min_value=1, max_value=int(1e7))
lengths = st.integers(min_value=1, max_value=int(1e4))

strands = st.sampled_from(VALID_GENOMIC_STRAND_INFO)

chromosomes = st.sampled_from([f"chr{e!s}" for e in list(range(1, 23)) + "X Y M".split()])


# dfs_zero_length_allowed = data_frames(
#     index=indexes(dtype=np.int64, min_size=0, unique=True, elements=starts),
#     columns=[
#         column("Chromosome", chromosomes),
#         column("Start", elements=starts),
#         column("End", elements=lengths),
#         column("Strand", strands),
#     ],
# )

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
        print("cmd " * 5)
        print(cmd)
        # ignoring the below line in bandit as only strings created by
        # the test suite is run here; no user input ever sought
        return subprocess.check_output(cmd, shell=True, executable="/bin/bash").decode()  # nosec  # nosec


@pytest.mark.bedtools()
@pytest.mark.parametrize("strand_behavior", ["same", "ignore"])
@settings(
    max_examples=max_examples,
    print_blob=True,
)
@given(gr=nonempty_pyranges(), gr2=nonempty_pyranges())  # pylint: disable=no-value-for-parameter
@reproduce_failure("6.92.1", b"AXicY2RABYxghA4AAHUABA==")
def test_set_intersect(gr, gr2, strand_behavior) -> None:
    print(max_examples)
    print(strand_behavior)
    set_intersect_command = "bedtools intersect {strand} -a <(sort -k1,1 -k2,2n {f1} | bedtools merge {strand} -c 4,5,6 -o first -i -) -b <(sort -k1,1 -k2,2n {f2} | bedtools merge {strand} -c 4,5,6 -o first -i -)"
    bedtools_result = run_bedtools(set_intersect_command, gr, gr2, strand_behavior)

    bedtools_df = read_bedtools_result_set_op(bedtools_result, strand_behavior=strand_behavior)

    result = gr.set_intersect(gr2, strand_behavior=strand_behavior)

    print("result:\n", result)
    print("bedtools_df\n", bedtools_df)
    pd.testing.assert_frame_equal(result, bedtools_df)
