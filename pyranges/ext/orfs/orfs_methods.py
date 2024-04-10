import numpy as np

import pyranges as pr
from pyranges.core.names import FRAME_COL, TEMP_CUMSUM_COL, TEMP_INDEX_COL, TEMP_LENGTH_COL
from pyranges.core.namespace_utils import decorate_to_pyranges_method
from pyranges.core.pyranges_helpers import mypy_ensure_pyranges
from pyranges.ext.orfs.extend_orfs import extend_orfs


def calculate_frame(self: pr.PyRanges, transcript_id: str | list[str], frame_col: str = "Frame") -> pr.PyRanges:
    """Calculate the frame of genomic intervals, assuming all are coding sequences (CDS), and add it as column.

    A stranded
    After this, the input Pyranges will contain an added "Frame" column, which determines the nucleotide of the CDS
    that is the first base of a codon.Resulting values are in range between 0 and 2 included.
    0 indicates that the first nucleotide of that interval is the first base of a codon,
    1 indicates the second base and 2 indicates the third base.
    While the 5'-most interval of each transcript has always 0 frame, the following ones may have any of these values.

    Parameters
    ----------
    self : PyRanges
        Input CDS intervals.

    transcript_id : str or list of str
        Column(s) to group by the intervals: coding exons belonging to the same transcript have the same values in this/these column(s).

    frame_col : str, default 'Frame'
        Name of the column to store the frame values.

    Returns
    -------
    PyRanges

    Examples
    --------
    >>> p = pr.PyRanges({"Chromosome": [1,1,1,2,2],
    ...                   "Strand": ["+","+","+","-","-"],
    ...                   "Start": [1,31,52,101,201],
    ...                   "End": [10,45,90,130,218],
    ...                   "transcript_id": ["t1","t1","t1","t2","t2"]})
    >>> p
      index  |      Chromosome  Strand      Start      End  transcript_id
      int64  |           int64  object      int64    int64  object
    -------  ---  ------------  --------  -------  -------  ---------------
          0  |               1  +               1       10  t1
          1  |               1  +              31       45  t1
          2  |               1  +              52       90  t1
          3  |               2  -             101      130  t2
          4  |               2  -             201      218  t2
    PyRanges with 5 rows, 5 columns, and 1 index columns.
    Contains 2 chromosomes and 2 strands.

    >>> p.orfs.calculate_frame(transcript_id=['transcript_id'])
      index  |      Chromosome  Strand      Start      End  transcript_id      Frame
      int64  |           int64  object      int64    int64  object             int64
    -------  ---  ------------  --------  -------  -------  ---------------  -------
          0  |               1  +               1       10  t1                     0
          1  |               1  +              31       45  t1                     0
          2  |               1  +              52       90  t1                     2
          3  |               2  -             101      130  t2                     2
          4  |               2  -             201      218  t2                     0
    PyRanges with 5 rows, 6 columns, and 1 index columns.
    Contains 2 chromosomes and 2 strands.

    """
    if not self.strand_valid:
        msg = "Strand must be valid to run calculate_frame."
        raise AssertionError(msg)

    gr = self.copy()

    # Column to save the initial index
    gr[TEMP_INDEX_COL] = np.arange(len(self))

    # Filtering for desired columns
    sorted_p = gr.get_with_loc_columns([TEMP_INDEX_COL, *self._by_to_list(transcript_id)])

    # Sorting by 5' (Intervals on + are sorted by ascending order and - are sorted by descending order)
    sorted_p = sorted_p.sort_by_5_prime_ascending_and_3_prime_descending()

    # Creating a column saving the length for the intervals (for selenoprofiles and ensembl)
    sorted_p[TEMP_LENGTH_COL] = sorted_p.lengths()

    # Creating a column saving the cumulative length for the intervals
    sorted_p[TEMP_CUMSUM_COL] = sorted_p.groupby(transcript_id)[TEMP_LENGTH_COL].cumsum()

    # Creating a frame column
    sorted_p[FRAME_COL] = (sorted_p[TEMP_CUMSUM_COL] - sorted_p[TEMP_LENGTH_COL]) % 3

    # Appending the Frame of sorted_p by the index of p
    sorted_p = mypy_ensure_pyranges(sorted_p.sort_values(by=TEMP_INDEX_COL))

    gr[frame_col] = sorted_p[FRAME_COL]

    return gr.drop_and_return(TEMP_INDEX_COL, axis=1)


class OrfsManager:
    """Namespace manager for statistical methods that accept PyRanges as first argument, accessed with pr.PyRanges.orfs.

    Additional methods are available at pyranges.orfs.
    """

    def __init__(self, p: pr.PyRanges) -> None:
        self.pyranges_instance = p

    calculate_frame = decorate_to_pyranges_method(calculate_frame)
    extend_orfs = decorate_to_pyranges_method(extend_orfs)
