Inspecting PyRanges
~~~~~~~~~~~~~~~~~~~

.. contents::
   :local:
   :depth: 2

String representation
---------------------

Print a PyRanges object for an overview of its data:

  >>> import pyranges as pr
  >>> gr = pr.example_data.chipseq
  >>> print(gr)
  index    |    Chromosome    Start      End        Name      Score    Strand
  int64    |    category      int64      int64      object    int64    category
  -------  ---  ------------  ---------  ---------  --------  -------  ----------
  0        |    chr8          28510032   28510057   U0        0        -
  1        |    chr7          107153363  107153388  U0        0        -
  2        |    chr5          135821802  135821827  U0        0        -
  3        |    chr14         19418999   19419024   U0        0        -
  ...      |    ...           ...        ...        ...       ...      ...
  16       |    chr9          120803448  120803473  U0        0        +
  17       |    chr6          89296757   89296782   U0        0        -
  18       |    chr1          194245558  194245583  U0        0        +
  19       |    chr8          57916061   57916086   U0        0        +
  PyRanges with 20 rows, 6 columns, and 1 index columns.
  Contains 15 chromosomes and 2 strands.

To obtain this representation, you can invoke the ``str`` builtin, e.g. with ``str(gr)``.

  >>> a = str(gr)
  >>> print(a)
  index    |    Chromosome    Start      End        Name      Score    Strand
  int64    |    category      int64      int64      object    int64    category
  -------  ---  ------------  ---------  ---------  --------  -------  ----------
  0        |    chr8          28510032   28510057   U0        0        -
  1        |    chr7          107153363  107153388  U0        0        -
  2        |    chr5          135821802  135821827  U0        0        -
  3        |    chr14         19418999   19419024   U0        0        -
  ...      |    ...           ...        ...        ...       ...      ...
  16       |    chr9          120803448  120803473  U0        0        +
  17       |    chr6          89296757   89296782   U0        0        -
  18       |    chr1          194245558  194245583  U0        0        +
  19       |    chr8          57916061   57916086   U0        0        +
  PyRanges with 20 rows, 6 columns, and 1 index columns.
  Contains 15 chromosomes and 2 strands.

Only a limited number of rows are displayed, which are taken from the top and bottom of the table.
You can change the number of rows displayed in any PyRanges using :func:`pyranges.options.set_options` as such:

  >>> pr.options.set_option('max_rows_to_show', 20)
  >>> gr
    index  |    Chromosome        Start        End  Name        Score  Strand
    int64  |    category          int64      int64  object      int64  category
  -------  ---  ------------  ---------  ---------  --------  -------  ----------
        0  |    chr8           28510032   28510057  U0              0  -
        1  |    chr7          107153363  107153388  U0              0  -
        2  |    chr5          135821802  135821827  U0              0  -
        3  |    chr14          19418999   19419024  U0              0  -
        4  |    chr12         106679761  106679786  U0              0  -
        5  |    chr21          40099618   40099643  U0              0  +
        6  |    chr8           22714402   22714427  U0              0  -
        7  |    chr19          19571102   19571127  U0              0  +
        8  |    chr3          140986358  140986383  U0              0  -
        9  |    chr10          35419784   35419809  U0              0  -
       10  |    chr4           98488749   98488774  U0              0  +
       11  |    chr11          22225193   22225218  U0              0  +
       12  |    chr1           38457520   38457545  U0              0  +
       13  |    chr1           80668132   80668157  U0              0  -
       14  |    chr2          152562484  152562509  U0              0  -
       15  |    chr4          153155301  153155326  U0              0  +
       16  |    chr9          120803448  120803473  U0              0  +
       17  |    chr6           89296757   89296782  U0              0  -
       18  |    chr1          194245558  194245583  U0              0  +
       19  |    chr8           57916061   57916086  U0              0  +
  PyRanges with 20 rows, 6 columns, and 1 index columns.
  Contains 15 chromosomes and 2 strands.

Let's reset display options to defaults:

  >>> pr.options.reset_options()



Detecting invalid PyRanges
--------------------------

The string representation of PyRanges shows useful information to detect data anomalies.

For example, intervals may have invalid lengths. Note that message at the bottom:

  >>> pr.PyRanges(dict(Chromosome='chr1', Start=[1, 10], End=[0, 20]))
    index  |    Chromosome      Start      End
    int64  |    object          int64    int64
  -------  ---  ------------  -------  -------
        0  |    chr1                1        0
        1  |    chr1               10       20
  PyRanges with 2 rows, 3 columns, and 1 index columns.
  Contains 1 chromosomes.
  Invalid ranges:
    * 1 intervals are empty or negative length (end <= start). See indexes: 0

Intervals may also be invalid because of NaN in their Start or End values:

  >>> pr.PyRanges(dict(Chromosome='chr1', Start=[None, 10], End=[0, 20]))
    index  |    Chromosome        Start      End
    int64  |    object          float64    int64
  -------  ---  ------------  ---------  -------
        0  |    chr1                nan        0
        1  |    chr1                 10       20
  PyRanges with 2 rows, 3 columns, and 1 index columns.
  Contains 1 chromosomes.
  Invalid ranges:
    * 1 starts or ends are nan. See indexes: 0

Or because they have negative Start/End values, see below. This can be remedied with
function :func:`clip_ranges <pyranges.PyRanges.clip_ranges>`.

  >>> pr.PyRanges(dict(Chromosome='chr1', Start=[1, -10], End=[11, 20]))
    index  |    Chromosome      Start      End
    int64  |    object          int64    int64
  -------  ---  ------------  -------  -------
        0  |    chr1                1       11
        1  |    chr1              -10       20
  PyRanges with 2 rows, 3 columns, and 1 index columns.
  Contains 1 chromosomes.
  Invalid ranges:
    * 1 starts or ends are < 0. See indexes: 1

A relatively common case is PyRanges objects that have a Strand column, but the strands are not valid genomic strands.
Note the warning in the last line of the string representation:

  >>> g = pr.PyRanges(dict(Chromosome='chr1', Start=[1, 1], End=[11, 20], Strand=['-', '#']))
  >>> g
    index  |    Chromosome      Start      End  Strand
    int64  |    object          int64    int64  object
  -------  ---  ------------  -------  -------  --------
        0  |    chr1                1       11  -
        1  |    chr1                1       20  #
  PyRanges with 2 rows, 4 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands (including non-genomic strands: #).

Non-valid strands can affect the functioning of many methods that have a ``use_strand`` parameter
(e.g. :func:`slice_ranges <pyranges.PyRanges.slice_ranges>`) or
a ``strand_behavior`` parameter (e.g. :func:`overlap <pyranges.PyRanges.overlap>`), because these parameters
by default are set to ``auto``, meaning that strand is considered only if it is valid.
Indeed, see that this subregion is calculated from the left limit, even for the interval on  the '-' strand:

  >>> g.slice_ranges(0, 3)
    index  |    Chromosome      Start      End  Strand
    int64  |    object          int64    int64  object
  -------  ---  ------------  -------  -------  --------
        0  |    chr1                1        4  -
        1  |    chr1                1        4  #
  PyRanges with 2 rows, 4 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands (including non-genomic strands: #).

When running the code above, you should get a warning message like this:

  .. code-block:: none

    UserWarning: slice_ranges: 'auto' use_strand treated as False due to invalid Strand values. Suppress this warning with use_strand=False
    g.slice_ranges(0, 3)

You can check whether a PyRanges object has valid Strand information with property
:func:`strand_valid <pyranges.PyRanges.strand_valid>`:

  >>> g.strand_valid
  False

To fix the invalid strands by turning them to '+',
use method :func:`make_strand_valid <pyranges.PyRanges.make_strand_valid>`:

  >>> g2 = g.make_strand_valid()
  >>> g2
    index  |    Chromosome      Start      End  Strand
    int64  |    object          int64    int64  object
  -------  ---  ------------  -------  -------  --------
        0  |    chr1                1       11  -
        1  |    chr1                1       20  +
  PyRanges with 2 rows, 4 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

Lastly, some operations may result in PyRanges with duplicated indices, which is shown in the
penultimate line of the string representation:

  >>> gr1= pr.PyRanges(dict(Chromosome='chr1', Start=[1], End=[100]))
  >>> gr2 = pr.PyRanges(dict(Chromosome='chr1', Start=[20, 50], End=[30, 60]))
  >>> gr3 = gr1.subtract_overlaps(gr2)
  >>> gr3
    index  |    Chromosome      Start      End
    int64  |    object          int64    int64
  -------  ---  ------------  -------  -------
        0  |    chr1                1       20
        0  |    chr1               30       50
        0  |    chr1               60      100
  PyRanges with 3 rows, 3 columns, and 1 index columns (with 2 index duplicates).
  Contains 1 chromosomes.

To remedy this, use pandas method ``reset_index``:

  >>> gr3 = gr3.reset_index(drop=True)
  >>> gr3
    index  |    Chromosome      Start      End
    int64  |    object          int64    int64
  -------  ---  ------------  -------  -------
        0  |    chr1                1       20
        1  |    chr1               30       50
        2  |    chr1               60      100
  PyRanges with 3 rows, 3 columns, and 1 index columns.
  Contains 1 chromosomes.

Column summary statistics
-------------------------
PyRanges columns are pandas Series, and they may be of different data types.
The types are shown in the header shown in their string representation (see above).
To see them all, use property ``dtypes`` like you do for dataframes:

  >>> gr.dtypes
  Chromosome    category
  Start            int64
  End              int64
  Name            object
  Score            int64
  Strand        category
  dtype: object

There are convenient methods inherited from pandas dataframes to inspect PyRanges objects, such as ``info``:

  >>> gr.info() # doctest: +NORMALIZE_WHITESPACE
  <class 'pyranges.core.pyranges_main.PyRanges'>
  RangeIndex: 20 entries, 0 to 19
  Data columns (total 6 columns):
   #   Column      Non-Null Count  Dtype
  ---  ------      --------------  -----
   0   Chromosome  20 non-null     category
   1   Start       20 non-null     int64
   2   End         20 non-null     int64
   3   Name        20 non-null     object
   4   Score       20 non-null     int64
   5   Strand      20 non-null     category
  dtypes: category(2), int64(3), object(1)
  memory usage: 1.6+ KB

On the other hand, ``describe`` reports aggregate metrics of numerical columns:

  >>> gr.describe()
                Start           End  Score
  count  2.000000e+01  2.000000e+01   20.0
  mean   8.320972e+07  8.320975e+07    0.0
  std    5.439939e+07  5.439939e+07    0.0
  min    1.941900e+07  1.941902e+07    0.0
  25%    3.369235e+07  3.369237e+07    0.0
  50%    8.498244e+07  8.498247e+07    0.0
  75%    1.245580e+08  1.245581e+08    0.0
  max    1.942456e+08  1.942456e+08    0.0
