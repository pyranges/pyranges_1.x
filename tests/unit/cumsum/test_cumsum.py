import pandas as pd

import pyranges as pr


def test_group_cumsum() -> None:
    """Test that group_cumsum works correctly with sorted/unsorted objects."""
    p = pr.PyRanges(
        {
            "Chromosome": [
                "chr2", "chr3", "chr2", "chr3", "chr3",
                "chr1", "chr2", "chr3", "chr3", "chr2",
                "chr2", "chr1", "chr1", "chr2", "chr3",
                "chr1", "chr2", "chr1", "chr3", "chr2"
            ],
            "Start": [
                785520, 89970, 626749, 168799, 372056,
                520458, 800809, 34720, 167104, 804868,
                59764, 972819, 235131, 4508, 866290,
                835796, 694284, 65536, 844945, 705072
            ],
            "End": [
                795577, 97475, 635187, 175999, 374410,
                521540, 807322, 41109, 167881, 805508,
                64323, 978796, 240207, 9362, 870401,
                845344, 701742, 73390, 851624, 712783
            ],
            "Strand": [
                "+", "-", "+", "+", "+",
                "+", "-", "+", "-", "-",
                "+", "+", "-", "+", "+",
                "-", "-", "-", "-", "+"
            ],
            "Label": [
                "grp0", "grp9", "grp2", "grp5", "grp5",
                "grp1", "grp7", "grp5", "grp9", "grp8",
                "grp2", "grp1", "grp3", "grp2", "grp5",
                "grp3", "grp7", "grp6", "grp9", "grp0"
            ]
        }    )

    def process_to_compare(df: pd.DataFrame) -> pd.DataFrame:
        return df[['Label', 'cs_end']].sort_values('cs_end')

    pcs=p.group_cumsum(group_by='Label', cumsum_end_column='cs_end', cumsum_start_column='cs_start')

    pcs_sorted = p.sort_ranges().group_cumsum(group_by='Label', cumsum_end_column='cs_end', cumsum_start_column='cs_start')

    pcs_sorted_lab = p.sort_ranges('Label').group_cumsum(group_by='Label', cumsum_end_column='cs_end',
                                              cumsum_start_column='cs_start')

    pd.testing.assert_frame_equal(process_to_compare(pcs), process_to_compare(pcs_sorted))
    pd.testing.assert_frame_equal(process_to_compare(pcs), process_to_compare(pcs_sorted_lab))



