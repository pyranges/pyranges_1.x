Operations on coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~

.. contents::
   :local:
   :depth: 2
   :caption: Contents


subseq, spl subseq,


Modifying interval coordinates
------------------------------
Interval coordinates can be directly modified like any Series in dataframes.

  >>> import pyranges as pr
  >>> e = pr.example_data.ensembl_gtf
  >>> e = e[e.Feature == "exon"].get_with_loc_columns('transcript_id')
  >>> e
    index  |      Chromosome    Start      End  Strand      transcript_id
    int64  |        category    int64    int64  category    object
  -------  ---  ------------  -------  -------  ----------  ---------------
        2  |               1    11868    12227  +           ENST00000456328
        3  |               1    12612    12721  +           ENST00000456328
        4  |               1    13220    14409  +           ENST00000456328
        5  |               1   112699   112804  -           ENST00000471248
        6  |               1   110952   111357  -           ENST00000471248
        8  |               1   133373   133723  -           ENST00000610542
        9  |               1   129054   129223  -           ENST00000610542
       10  |               1   120873   120932  -           ENST00000610542
  PyRanges with 8 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

We can modify a whole column at once:

  >>> e['Start'] += 5
  >>> e
    index  |      Chromosome    Start      End  Strand      transcript_id
    int64  |        category    int64    int64  category    object
  -------  ---  ------------  -------  -------  ----------  ---------------
        2  |               1    11873    12227  +           ENST00000456328
        3  |               1    12617    12721  +           ENST00000456328
        4  |               1    13225    14409  +           ENST00000456328
        5  |               1   112704   112804  -           ENST00000471248
        6  |               1   110957   111357  -           ENST00000471248
        8  |               1   133378   133723  -           ENST00000610542
        9  |               1   129059   129223  -           ENST00000610542
       10  |               1   120878   120932  -           ENST00000610542
  PyRanges with 8 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

Or we can modify a slice of the column:

  >>> e.loc[2:5, 'Start'] -= 5
  >>> e
    index  |      Chromosome    Start      End  Strand      transcript_id
    int64  |        category    int64    int64  category    object
  -------  ---  ------------  -------  -------  ----------  ---------------
        2  |               1    11868    12227  +           ENST00000456328
        3  |               1    12612    12721  +           ENST00000456328
        4  |               1    13220    14409  +           ENST00000456328
        5  |               1   112699   112804  -           ENST00000471248
        6  |               1   110957   111357  -           ENST00000471248
        8  |               1   133378   133723  -           ENST00000610542
        9  |               1   129059   129223  -           ENST00000610542
       10  |               1   120878   120932  -           ENST00000610542
  PyRanges with 8 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

Or use a boolean index:

  >>> e.loc[e.Strand == "+", "Start"] += 5
  >>> e
    index  |      Chromosome    Start      End  Strand      transcript_id
    int64  |        category    int64    int64  category    object
  -------  ---  ------------  -------  -------  ----------  ---------------
        2  |               1    11873    12227  +           ENST00000456328
        3  |               1    12617    12721  +           ENST00000456328
        4  |               1    13225    14409  +           ENST00000456328
        5  |               1   112699   112804  -           ENST00000471248
        6  |               1   110957   111357  -           ENST00000471248
        8  |               1   133378   133723  -           ENST00000610542
        9  |               1   129059   129223  -           ENST00000610542
       10  |               1   120878   120932  -           ENST00000610542
  PyRanges with 8 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

On the other hand, pyranges offer convenient and intuitive methods to modify coordinates, which deal
with the complexity of intervals and strands.

Subsequence operations
----------------------

Subsequence operations are operations that slice the intervals in a PyRanges object to obtain smaller intervals.
Intervals may be treated independently (default) or grouped in transcripts.

Method :func:`subsequence <pyranges.PyRanges.subsequence>` allows to
obtain subsequences by specifying the start and end position, in python notation:

  >>> #e.subsequence(0, 10)

Other slicing operations
------------------------

mention window, tile, tile_genome

