import pandas.testing
import pytest

import pyranges as pr


def test_split_chunks() -> None:
    df = pr.example_data.ensembl_gtf.get_with_loc_columns(["gene_id"])
    frame1, frame2 = pr.parallelism.split_df_into_chunks_without_splitting_groups(df, by=["gene_id"], nb_splits=2)
    assert len({*frame1["gene_id"]}) == 1
    assert len({*frame2["gene_id"]}) == 1
    assert len(frame1) + len(frame2) == len(df)


@pytest.mark.slow
def test_split_chunks_many() -> None:
    gr = pr.PyRanges(
        {"Chromosome": ["chr1"] * 1000, "Start": range(1000), "End": range(1, 1001), "GeneId": [*range(100)] * 10}
    )
    res = pr.parallelism.split_df_into_chunks_without_splitting_groups(gr, by=["GeneId"], nb_splits=37)
    assert sum(len(r) for r in res) == len(gr)
    seen_ids = set()
    for i, r in enumerate(res):
        new_ids = {*r["GeneId"].drop_duplicates()}
        assert new_ids.isdisjoint(seen_ids), i
        seen_ids.update(new_ids)


@pytest.mark.slow
def test_run_in_parallel() -> None:
    df = pr.example_data.ensembl_gtf.get_with_loc_columns(["gene_id"])
    dfs = [gdf for _, gdf in df.groupby("gene_id")]
    res = pr.parallelism.run_in_parallel(
        function=lambda x: x.merge_overlaps(by="gene_id", count_col="Counts"),
        dfs=dfs,
        nb_cpu=2,
    ).reset_index(drop=True)
    expected_result = df.merge_overlaps(by="gene_id", count_col="Counts")
    pandas.testing.assert_frame_equal(res, expected_result)
