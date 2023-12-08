import logging
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import pandas as pd
from pandas.core.frame import DataFrame
from pandas.core.series import Series

from pyranges.names import END_COL, FORWARD_STRAND, START_COL

if TYPE_CHECKING:
    import pyfaidx  # type: ignore[import]

    from pyranges import PyRanges


logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)
LOGGER.Formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s %(message)s")
LOGGER.setLevel(logging.INFO)

def get_sequence(
    gr: "PyRanges",
    path: Path | None = None,
    pyfaidx_fasta: Optional["pyfaidx.Fasta"] = None,
) -> Series:
    r"""Get the sequence of the intervals from a fasta file.

    Parameters
    ----------
    gr : PyRanges
        Coordinates.

    path : Path
        Path to fasta file. It will be indexed using pyfaidx if an index is not found

    pyfaidx_fasta : pyfaidx.Fasta
        Alternative method to provide fasta target, as a pyfaidx.Fasta object


    Returns
    -------
    Series

        Sequences, one per interval.

    Note
    ----

    This function requires the library pyfaidx, it can be installed with
    ``conda install -c bioconda pyfaidx`` or ``pip install pyfaidx``.

    Sorting the PyRanges is likely to improve the speed.
    Intervals on the negative strand will be reverse complemented.

    Warning
    -------

    Note that the names in the fasta header and gr must be the same.

    See Also
    --------
    get_transcript_sequence : obtain mRNA sequences, by joining exons belonging to the same transcript


    Examples
    --------
    >>> import pyranges as pr
    >>> gr = pr.PyRanges({"Chromosome": ["chr1", "chr1"],
    ...                   "Start": [5, 0], "End": [8, 5],
    ...                   "Strand": ["+", "-"]})

    >>> gr
      index  |    Chromosome      Start      End  Strand
      int64  |    object          int64    int64  object
    -------  ---  ------------  -------  -------  --------
          0  |    chr1                5        8  +
          1  |    chr1                0        5  -
    PyRanges with 2 rows, 4 columns, and 1 index columns.
    Contains 1 chromosomes and 2 strands.

    >>> tmp_handle = open("temp.fasta", "w+")
    >>> _ = tmp_handle.write(">chr1\n")
    >>> _ = tmp_handle.write("GTAATCAT\n")
    >>> tmp_handle.close()

    >>> seq = pr.get_sequence(gr, "temp.fasta")

    >>> seq
    0      CAT
    1    ATTAC
    dtype: object

    >>> gr.col.seq = seq
    >>> gr
      index  |    Chromosome      Start      End  Strand    seq
      int64  |    object          int64    int64  object    object
    -------  ---  ------------  -------  -------  --------  --------
          0  |    chr1                5        8  +         CAT
          1  |    chr1                0        5  -         ATTAC
    PyRanges with 2 rows, 5 columns, and 1 index columns.
    Contains 1 chromosomes and 2 strands.
    """
    try:
        import pyfaidx  # type: ignore[import]
    except ImportError:
        LOGGER.exception(
            "pyfaidx must be installed to get fasta sequences. Use `conda install -c bioconda pyfaidx` or `pip install pyfaidx` to install it.",
        )
        sys.exit(1)

    if pyfaidx_fasta is None:
        if path is None:
            msg = "ERROR get_sequence : you must provide a fasta path or pyfaidx_fasta object"
            raise ValueError(msg)
        pyfaidx_fasta = pyfaidx.Fasta(path, read_ahead=int(1e5))

    seqs = []
    use_strand = gr.strand_values_valid
    for key, df in gr.groupby(gr.location_cols_include_strand_only_if_valid):
        chromosome, strand = key + (() if use_strand else (FORWARD_STRAND,))
        _fasta = pyfaidx_fasta[chromosome]
        forward_strand = strand == FORWARD_STRAND
        for start, end in zip(df[START_COL], df[END_COL], strict=True):
            seq = _fasta[start:end]
            seqs.append(seq.seq if forward_strand else (-seq).seq)

    return pd.concat([pd.Series(s) for s in seqs]).reset_index(drop=True).squeeze()


def get_transcript_sequence(
    gr: "PyRanges",
    group_by: str,
    path: Path | None = None,
    pyfaidx_fasta: Optional["pyfaidx.Fasta"] = None,
) -> DataFrame:
    r"""Get the sequence of mRNAs, e.g. joining intervals corresponding to exons of the same transcript.

    Parameters
    ----------
    gr : PyRanges
        Coordinates.

    group_by : str or list of str
        intervals are grouped by this/these ID column(s): these are exons belonging to same transcript

    path : Optional Path
        Path to fasta file. It will be indexed using pyfaidx if an index is not found

    pyfaidx_fasta : pyfaidx.Fasta
        Alternative method to provide fasta target, as a pyfaidx.Fasta object


    Returns
    -------
    DataFrame

        Pandas DataFrame with a column for Sequence, plus ID column(s) provided with "group_by"

    Note
    ----

    This function requires the library pyfaidx, it can be installed with
    ``conda install -c bioconda pyfaidx`` or ``pip install pyfaidx``.

    Sorting the PyRanges is likely to improve the speed.
    Intervals on the negative strand will be reverse complemented.

    Warning
    -------

    Note that the names in the fasta header and gr must be the same.

    See Also
    --------
    get_sequence : obtain sequence of single intervals


    Examples
    --------
    >>> import pyranges as pr
    >>> gr = pr.PyRanges({"Chromosome": ['chr1', 'chr1', 'chr1'],
    ...                   "Start": [0, 9, 18], "End": [4, 13, 21],
    ...                   "Strand":['+', '-', '-'],
    ...                   "transcript": ['t1', 't2', 't2']})
    >>> gr
      index  |    Chromosome      Start      End  Strand    transcript
      int64  |    object          int64    int64  object    object
    -------  ---  ------------  -------  -------  --------  ------------
          0  |    chr1                0        4  +         t1
          1  |    chr1                9       13  -         t2
          2  |    chr1               18       21  -         t2
    PyRanges with 3 rows, 5 columns, and 1 index columns.
    Contains 1 chromosomes and 2 strands.

    >>> tmp_handle = open("temp.fasta", "w+")
    >>> _ = tmp_handle.write(">chr1\n")
    >>> _ = tmp_handle.write("AAACCCTTTGGGAAACCCTTTGGG\n")
    >>> tmp_handle.close()

    >>> seq = pr.get_transcript_sequence(gr, path="temp.fasta", group_by='transcript')
    >>> seq
      transcript Sequence
    0         t1     AAAC
    1         t2  AAATCCC

    To write to a file in fasta format:
    >>> with open('outfile.fasta', 'w') as fw:
    ...     nchars=60
    ...     for row in seq.itertuples():
    ...         s = '\\n'.join([ row.Sequence[i:i+nchars] for i in range(0, len(row.Sequence), nchars)])
    ...         _bytes_written = fw.write(f'>{row.transcript}\\n{s}\\n')
    """
    gr = gr.sort_by_5_prime_ascending_and_3_prime_descending() if gr.strand_values_valid else gr.sort_by_position()

    gr.col["Sequence"] = get_sequence(gr, path=path, pyfaidx_fasta=pyfaidx_fasta)

    return gr.groupby(group_by, as_index=False).agg({"Sequence": "".join})
