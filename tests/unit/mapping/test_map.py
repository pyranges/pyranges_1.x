import pandas as pd
import pyranges as pr

def test_map_to_global() -> None:
    """Test that map_to_global works correctly with sorted/unsorted objects."""

    def process_to_compare(df: pd.DataFrame) -> pd.DataFrame:
        return df.sort_ranges()

    rh = pr.example_data.rfam_hits
    rh = rh[['target_name', 'seq_from', 'seq_to', 'strand', 'query_name', 'mdl_from', 'mdl_to']]
    rh = rh.rename(columns={'target_name': 'Chromosome', 'seq_from': 'Start', 'seq_to': 'End', 'strand': 'Strand'})
    rh = pr.PyRanges(rh)
    rh['Start'] -= 1  # convert to 0-based coordinates
    rh.loc[rh.Strand == '-', ['Start', 'End']] = rh.loc[rh.Strand == '-', ['End', 'Start']].values
    gr = pr.example_data.ncbi_gff
    gre = (gr[gr.Feature == 'exon']).get_with_loc_columns('Parent')
    rh = rh.remove_nonloc_columns()

    rhg_def = rh.map_to_global(gre, global_on='Parent', keep_id=True, keep_loc=True)
    rhg_sl1 = rh.sort_ranges().map_to_global(gre, global_on='Parent', keep_id=True, keep_loc=True)
    rhg_sg1 = rh.map_to_global(gre.sort_ranges(), global_on='Parent', keep_id=True, keep_loc=True)
    rhg_sg2 = rh.map_to_global(gre.sort_ranges('Parent'), global_on='Parent', keep_id=True, keep_loc=True)

    assert process_to_compare(rhg_def).equals(process_to_compare(rhg_sl1))
    assert process_to_compare(rhg_def).equals(process_to_compare(rhg_sg1))
    assert process_to_compare(rhg_def).equals(process_to_compare(rhg_sg2))

    # gre data is bad!