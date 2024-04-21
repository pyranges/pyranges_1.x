How-to pages
============

Loading/Creating PyRanges
~~~~~~~~~~~~~~~~~~~~~~~~~


A PyRanges object can be built like a Pandas DataFrame, but genomic location columns (Chromosome, Start, End) are
mandatory. Refer to https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html for more information on how
to create a DataFrame. In alternative, PyRanges can be read from a file in bed, gtf, gff3 or bam format.

PyRanges are created in the following ways:

#. from a pandas dataframe
#. from a dictionary with the column names as keys and iterables as values
#. from a file, using pyranges readers
#. concatenating existing PyRanges objects

From a DataFrame
----------------

If you instantiate a PyRanges object from a dataframe, it should at least contain the columns Chromosome, Start and End.
Coordinates follow the python standard (0-based, start included, end excluded). A column called Strand is optional.
Any other columns in the dataframe are carried over as metadata.

  >>> import pandas as pd, pyranges as pr, numpy as np
  >>> df=pd.DataFrame(
  ... {'Chromosome':['chr1', 'chr1', 'chr1', 'chr3'],
  ...  'Start': [5, 20, 80, 10],
  ...  'End':   [10, 28, 95, 38],
  ...  'Strand':['+', '+', '-', '+'],
  ...  'title': ['a', 'b', 'c', 'd']}
  ... )
  >>> df
    Chromosome  Start  End Strand title
  0       chr1      5   10      +     a
  1       chr1     20   28      +     b
  2       chr1     80   95      -     c
  3       chr3     10   38      +     d


To instantiate PyRanges from a dataframe, provide it as argument to the PyRanges constructor:

  >>> p=pr.PyRanges(df)
  >>> p
    index  |    Chromosome      Start      End  Strand    title
    int64  |    object          int64    int64  object    object
  -------  ---  ------------  -------  -------  --------  --------
        0  |    chr1                5       10  +         a
        1  |    chr1               20       28  +         b
        2  |    chr1               80       95  -         c
        3  |    chr3               10       38  +         d
  PyRanges with 4 rows, 5 columns, and 1 index columns.
  Contains 2 chromosomes and 2 strands.


From a Dictionary
-----------------

You can instantiate a PyRanges object using dictionary:

  >>> gr = pr.PyRanges({'Chromosome': ['chr1', 'chr1', 'chr1', 'chr3'],
  ...                   'Start': [5, 20, 80, 10],
  ...                   'End': [10, 28, 95, 38],
  ...                   'Strand': ['+', '+', '-', '+'],
  ...                   'title': ['a', 'b', 'c', 'd']})
  >>> gr
    index  |    Chromosome      Start      End  Strand    title
    int64  |    object          int64    int64  object    object
  -------  ---  ------------  -------  -------  --------  --------
        0  |    chr1                5       10  +         a
        1  |    chr1               20       28  +         b
        2  |    chr1               80       95  -         c
        3  |    chr3               10       38  +         d
  PyRanges with 4 rows, 5 columns, and 1 index columns.
  Contains 2 chromosomes and 2 strands.

Both the ``{}`` and the ``dict`` constructors can be used to create a dictionary:

  >>> gr2 = pr.PyRanges(dict(Chromosome=['chr10', 'chr10', 'chr1', 'chr3'],
  ...                       Start=[55, 250, 80, 100],
  ...                       End=[150, 258, 95, 380]))
  >>> gr2
    index  |    Chromosome      Start      End
    int64  |    object          int64    int64
  -------  ---  ------------  -------  -------
        0  |    chr10              55      150
        1  |    chr10             250      258
        2  |    chr1               80       95
        3  |    chr3              100      380
  PyRanges with 4 rows, 3 columns, and 1 index columns.
  Contains 3 chromosomes.


As in the creation of a Dataframe, each list in the dictionary must have the same length.
Also, it may be any iterable, including pd.Series or np.array.
If a string or scalar is provided, it is broadcasted to the length of the other columns:

  >>> gr3 = pr.PyRanges(dict(Chromosome=pd.Series(['chr10', 'chr10', 'chr1', 'chr3']),
  ...                       Start=np.array([55, 250, 80, 100]),
  ...                       End=[150, 258, 95, 380],
  ...                       Strand='+'))
  >>> gr3
    index  |    Chromosome      Start      End  Strand
    int64  |    object          int64    int64  object
  -------  ---  ------------  -------  -------  --------
        0  |    chr10              55      150  +
        1  |    chr10             250      258  +
        2  |    chr1               80       95  +
        3  |    chr3              100      380  +
  PyRanges with 4 rows, 4 columns, and 1 index columns.
  Contains 3 chromosomes and 1 strands.


Loading from a file
-------------------

The pyranges library can create PyRanges from gff3 common file formats, namely gtf/gff, gff3, bed and bam (see
:func:`read_bed <pyranges.read_bed>`, :func:`read_gtf <pyranges.read_gtf>`,
:func:`read_gff3 <pyranges.read_gff3>`, :func:`read_bam <pyranges.read_bam>`).
The documentation of readers is available in the :doc:`pyranges module <pyranges_module>`.
Note that these files may encode coordinates with different conventions (e.g. GTF: 1-based, start and end included).
When instancing a PyRanges object they are converted to the python convention.

  >>> ensembl_path = pr.example_data.files['ensembl.gtf']  # example file
  >>> gr = pr.read_gtf(ensembl_path)
  >>> gr
  index    |    Chromosome    Source    Feature     Start    End      Score     Strand      Frame     gene_id          ...
  int64    |    category      object    category    int64    int64    object    category    object    object           ...
  -------  ---  ------------  --------  ----------  -------  -------  --------  ----------  --------  ---------------  -----
  0        |    1             havana    gene        11868    14409    .         +           .         ENSG00000223972  ...
  1        |    1             havana    transcript  11868    14409    .         +           .         ENSG00000223972  ...
  2        |    1             havana    exon        11868    12227    .         +           .         ENSG00000223972  ...
  3        |    1             havana    exon        12612    12721    .         +           .         ENSG00000223972  ...
  ...      |    ...           ...       ...         ...      ...      ...       ...         ...       ...              ...
  8        |    1             ensembl   transcript  120724   133723   .         -           .         ENSG00000238009  ...
  9        |    1             ensembl   exon        133373   133723   .         -           .         ENSG00000238009  ...
  10       |    1             ensembl   exon        129054   129223   .         -           .         ENSG00000238009  ...
  11       |    1             ensembl   exon        120873   120932   .         -           .         ENSG00000238009  ...
  PyRanges with 12 rows, 23 columns, and 1 index columns. (14 columns not shown: "gene_version", "gene_name", "gene_source", ...).
  Contains 1 chromosomes and 2 strands.


To read bam files the optional bamread-library must be installed. Use::

    conda install -c bioconda bamread

or::

    pip install bamread

to install it.

  >>> bam_path = pr.example_data.files['smaller.bam']  # example file
  >>> gr4 = pr.read_bam(bam_path)
  >>> gr4
  index    |    Chromosome    Start     End       Strand      Flag
  int64    |    category      int64     int64     category    uint16
  -------  ---  ------------  --------  --------  ----------  --------
  0        |    chr1          887771    887796    -           16
  1        |    chr1          994660    994685    -           16
  2        |    chr1          1041102   1041127   +           0
  3        |    chr1          1770383   1770408   -           16
  ...      |    ...           ...       ...       ...         ...
  96       |    chr1          18800901  18800926  +           0
  97       |    chr1          18800901  18800926  +           0
  98       |    chr1          18855123  18855148  -           16
  99       |    chr1          19373470  19373495  +           0
  PyRanges with 100 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

read_bam takes the arguments ``sparse``, ``mapq``, ``required_flag``, ``filter_flag``,
which have the default values True, 0, 0 and 1540, respectively.
With sparse True, only the columns ``['Chromosome', 'Start', 'End', 'Strand', 'Flag']``
are fetched. Setting sparse to False additionally gives you the columns
``['QueryStart', 'QueryEnd', 'QuerySequence', 'Name', 'Cigar', 'Quality']``, but is more time and memory-consuming:

  >>> pr.read_bam(bam_path, sparse=False)
  index    |    Chromosome    Start     End       Strand      Flag      QueryStart    QueryEnd    QuerySequence    ...
  int64    |    category      int64     int64     category    uint16    int64         int64       object           ...
  -------  ---  ------------  --------  --------  ----------  --------  ------------  ----------  ---------------  -----
  0        |    chr1          887771    887796    -           16        0             25          None             ...
  1        |    chr1          994660    994685    -           16        0             25          None             ...
  2        |    chr1          1041102   1041127   +           0         0             25          None             ...
  3        |    chr1          1770383   1770408   -           16        0             25          None             ...
  ...      |    ...           ...       ...       ...         ...       ...           ...         ...              ...
  96       |    chr1          18800901  18800926  +           0         0             25          None             ...
  97       |    chr1          18800901  18800926  +           0         0             25          None             ...
  98       |    chr1          18855123  18855148  -           16        0             25          None             ...
  99       |    chr1          19373470  19373495  +           0         0             25          None             ...
  PyRanges with 100 rows, 11 columns, and 1 index columns. (3 columns not shown: "Name", "Cigar", "Quality").
  Contains 1 chromosomes and 2 strands.

To load tabular file in any format, you can use pandas ``read_csv`` method and then pass the resulting dataframe to the
PyRanges constructor. Be aware of the coordinate convention of the file you load, and make sure that the
dataframe has aptly named columns.

From the concatenation of PyRanges
----------------------------------

Analogously to ``pandas.concat``, :func:`pyranges.concat` can be used to concatenate PyRanges objects, i.e.
stack rows of two or more PyRanges to create a new PyRanges object.

  >>> gr1 = pr.PyRanges({'Chromosome': ['chr1', 'chr1', 'chr1', 'chr3'],
  ...                    'Start': [5, 20, 80, 10],
  ...                    'End': [10, 28, 95, 38],
  ...                    'Strand': ['+', '+', '-', '+'],
  ...                    'title': ['a', 'b', 'c', 'd']})
  >>> gr2 = pr.PyRanges({'Chromosome': ['chr1', 'chr1', 'chr1', 'chr3'],
  ...                    'Start': [5, 20, 80, 10],
  ...                    'End': [10, 28, 95, 38],
  ...                    'Strand': ['+', '+', '-', '+'],
  ...                    'title': ['a', 'b', 'c', 'd']})
  >>> gr3 = pr.concat([gr1, gr2])
  >>> gr3
    index  |    Chromosome      Start      End  Strand    title
    int64  |    object          int64    int64  object    object
  -------  ---  ------------  -------  -------  --------  --------
        0  |    chr1                5       10  +         a
        1  |    chr1               20       28  +         b
        2  |    chr1               80       95  -         c
        3  |    chr3               10       38  +         d
        0  |    chr1                5       10  +         a
        1  |    chr1               20       28  +         b
        2  |    chr1               80       95  -         c
        3  |    chr3               10       38  +         d
  PyRanges with 8 rows, 5 columns, and 1 index columns (with 4 index duplicates).
  Contains 2 chromosomes and 2 strands.

Note that this may result in index duplicates, which can be remedied by pandas ``reset_index`` method.

  >>> pr.concat([gr1, gr2]).reset_index(drop=True)
    index  |    Chromosome      Start      End  Strand    title
    int64  |    object          int64    int64  object    object
  -------  ---  ------------  -------  -------  --------  --------
        0  |    chr1                5       10  +         a
        1  |    chr1               20       28  +         b
        2  |    chr1               80       95  -         c
        3  |    chr3               10       38  +         d
        4  |    chr1                5       10  +         a
        5  |    chr1               20       28  +         b
        6  |    chr1               80       95  -         c
        7  |    chr3               10       38  +         d
  PyRanges with 8 rows, 5 columns, and 1 index columns.
  Contains 2 chromosomes and 2 strands.

Obtaining PyRanges for testing purposes
----------------------------------------

For testing purposes, pyranges provides some data in ``pr.example_data``. See an overview with:

  >>> pr.example_data
  Available example data:
  -----------------------
  example_data.chipseq            : Example ChIP-seq data.
  example_data.chipseq_background : Example ChIP-seq data.
  example_data.chromsizes         : Example chromsizes data (hg19).
  example_data.ensembl_gtf        : Example gtf file from Ensembl.
  example_data.f1                 : Example bed file.
  example_data.f2                 : Example bed file.
  example_data.aorta              : Example ChIP-seq data.
  example_data.aorta2             : Example ChIP-seq data.
  example_data.ncbi_gff           : Example NCBI GFF data.
  example_data.ncbi_fasta         : Example NCBI fasta.

You can load the data with this syntax:

  >>> cs = pr.example_data.chipseq
  >>> cs
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

On the on other hand, you can create random intervals using :func:`pyranges.random`.
By default, the data refers to the human genome (hg19):

  >>> pr.random(n=5, length=50, seed=123)
    index  |    Chromosome        Start        End  Strand
    int64  |    object            int64      int64  object
  -------  ---  ------------  ---------  ---------  --------
        0  |    chr12          46554011   46554061  +
        1  |    chr1          202415019  202415069  -
        2  |    chr3           89385747   89385797  -
        3  |    chr3          182842974  182843024  -
        4  |    chr3           89004288   89004338  -
  PyRanges with 5 rows, 4 columns, and 1 index columns.
  Contains 3 chromosomes and 2 strands.


Writing to disk
~~~~~~~~~~~~~~~


The PyRanges can be written to several formats, namely csv, gtf, gff3 and bigwig.
If no path-argument is given, the string representation of the data is returned. (It may potentially be very large.)
If a path is given, it is taken as the path to the file to be written; in this case, the return value is the object
itself, to allow inserting write methods into method call chains.


Writing in tabular formats
--------------------------

Tabular formats such as csv, gtf, gff3 are the most popular for genomic annotations.
You can readily write them using the correspondent methods (see
:func:`to_bed <pyranges.PyRanges.to_bed>`,
:func:`to_gtf <pyranges.PyRanges.to_gtf>`,
:func:`to_gff3 <pyranges.PyRanges.to_gff3>`).

  >>> import pyranges as pr
  >>> gr = pr.example_data.chipseq
  >>> gr.to_gtf("chipseq.gtf")
  >>> #file chipseq.gtf has been created



Methods to_gff3 and to_gtf have a default mapping of PyRanges columns to GFF/GTF fields.
All extra ("metadata") columns are put in the last field:

  >>> gr['Label']='something' # bug below; to be updated with to_gtf is fixed
  >>> print(gr.head().to_gtf()) # doctest: +NORMALIZE_WHITESPACE
  chr8	.	.	28510033	28510057	0	-	.	Name=U0Label=something
  chr7	.	.	107153364	107153388	0	-	.	Name=U0Label=something
  chr5	.	.	135821803	135821827	0	-	.	Name=U0Label=something
  chr14	.	.	19419000	19419024	0	-	.	Name=U0Label=something
  chr12	.	.	106679762	106679786	0	-	.	Name=U0Label=something

Such mapping, as well as which attribute(s) are included as last field, can be altered. See the API for details.

The csv format is the most flexible, as it allows for any column to be included, and any separator to be used.
The method ``to_csv`` is directly inherited by pandas, so search for its API for details.
Remember that ``to_csv`` will not alter coordinates, so the output
will have the same pythonic convention as PyRanges:

The ``to_csv`` method takes the arguments header and sep:

  >>> print(gr.drop(['Label'], axis=1).head().to_csv(sep="\t", header=False, index=False)) # doctest: +NORMALIZE_WHITESPACE
  chr8	28510032	28510057	U0	0	-
  chr7	107153363	107153388	U0	0	-
  chr5	135821802	135821827	U0	0	-
  chr14	19418999	19419024	U0	0	-
  chr12	106679761	106679786	U0	0	-
  <BLANKLINE>


The `bigwig <http://genome.ucsc.edu/goldenPath/help/bigWig.html>`_ format differs substantially from
the formats above. Bigwig is a binary format, and it is typically used for large continuous quantitative
data along a genome sequence.

The pyranges library can also create bigwigs, but it needs the library pybigwig which is not installed by default.

Use this to install it::

	pip install pybigwig

The bigwig writer needs to know the chromosome sizes, e.g. provided as a dictionary {chromosome_name: size}.
You can also derive chromosome sizes from a fasta file using pyfaidx (see above to install it).

.. doctest::

  >>> import pyfaidx
  >>> chromsizes=pyfaidx.Fasta('your_genome.fa') # doctest: +SKIP


Once you obtained chromosome sizes, you are ready to write your PyRanges object to a bigwig file:

  >>> gr.to_bigwig("chipseq.bw", chromsizes) # doctest: +SKIP
  >>> # file chipseq.bw has been created

Bigwig is typically used to represent a coverage of some type.
To compute it from an arbitrary value column, use the value_col argument. See the API for additional options.
If you want to write one bigwig for each strand, you need to do it manually.


  >>> gr.loci["+"].to_bigwig("chipseq_plus.bw", chromsizes) # doctest: +SKIP
  >>> gr.loci["-"].to_bigwig("chipseq_minus.bw", chromsizes) # doctest: +SKIP


Inspecting PyRanges
~~~~~~~~~~~~~~~~~~~

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

Other convenient pandas methods are available to inspect PyRanges objects, such as:

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

Accessing data
~~~~~~~~~~~~~~

Selecting rows
--------------

Indexing with iloc, loc
+++++++++++++++++++++++

PyRanges inherits all the indexing and slicing capabilities of pandas, e.g. boolean Series indexing,
``iloc``, ``loc``, ``at``, ``iat``.
Note that these methods return a view, not a copy, with the caveats that it implies.
See the pandas documentation for details.
Briefly, to avoid ambiguity it is best to explicitly call ``copy`` if you want an object to not be linked
to the original object from which it was extracted. For example:

  >>> gr = pr.example_data.aorta
  >>> gr
  index    |    Chromosome    Start    End      Name      Score    Strand
  int64    |    category      int64    int64    object    int64    category
  -------  ---  ------------  -------  -------  --------  -------  ----------
  0        |    chr1          9916     10115    H3K27me3  5        -
  1        |    chr1          9939     10138    H3K27me3  7        +
  2        |    chr1          9951     10150    H3K27me3  8        -
  3        |    chr1          9953     10152    H3K27me3  5        +
  ...      |    ...           ...      ...      ...       ...      ...
  7        |    chr1          10127    10326    H3K27me3  1        -
  8        |    chr1          10241    10440    H3K27me3  6        -
  9        |    chr1          10246    10445    H3K27me3  4        +
  10       |    chr1          110246   110445   H3K27me3  1        +
  PyRanges with 11 rows, 6 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

  >>> sgr = gr.iloc[0:3].copy()
  >>> sgr
    index  |    Chromosome      Start      End  Name        Score  Strand
    int64  |    category        int64    int64  object      int64  category
  -------  ---  ------------  -------  -------  --------  -------  ----------
        0  |    chr1             9916    10115  H3K27me3        5  -
        1  |    chr1             9939    10138  H3K27me3        7  +
        2  |    chr1             9951    10150  H3K27me3        8  -
  PyRanges with 3 rows, 6 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

  >>> sgr['Score'] = 100  # does not modify gr
  >>> gr.head(3)
    index  |    Chromosome      Start      End  Name        Score  Strand
    int64  |    category        int64    int64  object      int64  category
  -------  ---  ------------  -------  -------  --------  -------  ----------
        0  |    chr1             9916    10115  H3K27me3        5  -
        1  |    chr1             9939    10138  H3K27me3        7  +
        2  |    chr1             9951    10150  H3K27me3        8  -
  PyRanges with 3 rows, 6 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

On the other hand, to modify even a few lines of a PyRanges object,
use ``loc`` (label-based indexing) or ``iloc`` (positional indexing) on the whole object, not on a view:

  >>> gr.loc[0:2, 'Score'] = 100  # modifies gr
  >>> gr.head(5)
    index  |    Chromosome      Start      End  Name        Score  Strand
    int64  |    category        int64    int64  object      int64  category
  -------  ---  ------------  -------  -------  --------  -------  ----------
        0  |    chr1             9916    10115  H3K27me3      100  -
        1  |    chr1             9939    10138  H3K27me3      100  +
        2  |    chr1             9951    10150  H3K27me3      100  -
        3  |    chr1             9953    10152  H3K27me3        5  +
        4  |    chr1             9978    10177  H3K27me3        7  -
  PyRanges with 5 rows, 6 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

Using boolean indexers
++++++++++++++++++++++
Analogous example using a boolean indexer, already seen in the tutorial:

  >>> gr.loc[gr['Score'] < 6, 'Score'] = -10  # modifies gr
  >>> gr.head(5)
    index  |    Chromosome      Start      End  Name        Score  Strand
    int64  |    category        int64    int64  object      int64  category
  -------  ---  ------------  -------  -------  --------  -------  ----------
        0  |    chr1             9916    10115  H3K27me3      100  -
        1  |    chr1             9939    10138  H3K27me3      100  +
        2  |    chr1             9951    10150  H3K27me3      100  -
        3  |    chr1             9953    10152  H3K27me3      -10  +
        4  |    chr1             9978    10177  H3K27me3        7  -
  PyRanges with 5 rows, 6 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.


In pandas, these logical operators can be employed with boolean Series:

* "&" =  element-wise AND operator
* "|" = element-wise OR operator
* "~" = NOT operator, inverts the values of the Series on its right

When using logical operators, make sure to parenthesize properly.

Let's get the + intervals with Score<8 starting before 10,000 or ending after 100,000:

  >>> gr[ (gr.Score<8) & (gr.Strand=='+') &
  ...     ((gr.Start<10000) | (gr.End>100000)) ]
    index  |    Chromosome      Start      End  Name        Score  Strand
    int64  |    category        int64    int64  object      int64  category
  -------  ---  ------------  -------  -------  --------  -------  ----------
        3  |    chr1             9953    10152  H3K27me3      -10  +
       10  |    chr1           110246   110445  H3K27me3      -10  +
  PyRanges with 2 rows, 6 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.

Let's invert the selection:
  >>> gr[~(
  ...      (gr.Score<8) & (gr.Strand=='+') &
  ...      ((gr.Start<10000) | (gr.End>100000)) )]
  index    |    Chromosome    Start    End      Name      Score    Strand
  int64    |    category      int64    int64    object    int64    category
  -------  ---  ------------  -------  -------  --------  -------  ----------
  0        |    chr1          9916     10115    H3K27me3  100      -
  1        |    chr1          9939     10138    H3K27me3  100      +
  2        |    chr1          9951     10150    H3K27me3  100      -
  4        |    chr1          9978     10177    H3K27me3  7        -
  ...      |    ...           ...      ...      ...       ...      ...
  6        |    chr1          10024    10223    H3K27me3  -10      +
  7        |    chr1          10127    10326    H3K27me3  -10      -
  8        |    chr1          10241    10440    H3K27me3  6        -
  9        |    chr1          10246    10445    H3K27me3  -10      +
  PyRanges with 9 rows, 6 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

Using PyRanges .loci
++++++++++++++++++++

As seen in the tutorial, PyRanges also provides the method :func:`loci <pyranges.PyRanges.loci>`
to select rows by genomic region:

  >>> gr2 = pr.example_data.aorta2.sort_ranges()
  >>> gr2
  index    |    Chromosome    Start    End      Name      Score    Strand
  int64    |    category      int64    int64    object    int64    category
  -------  ---  ------------  -------  -------  --------  -------  ----------
  1        |    chr1          10073    10272    Input     1        +
  5        |    chr1          10280    10479    Input     1        +
  6        |    chr1          16056    16255    Input     1        +
  7        |    chr1          16064    16263    Input     1        +
  ...      |    ...           ...      ...      ...       ...      ...
  4        |    chr1          10149    10348    Input     1        -
  3        |    chr1          10082    10281    Input     1        -
  2        |    chr1          10079    10278    Input     1        -
  0        |    chr1          9988     10187    Input     1        -
  PyRanges with 10 rows, 6 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

Various syntaxes are accepted, see its API. For example:

  >>> gr2.loci['-'] # get all rows with strand '-'
    index  |    Chromosome      Start      End  Name        Score  Strand
    int64  |    category        int64    int64  object      int64  category
  -------  ---  ------------  -------  -------  --------  -------  ----------
        9  |    chr1            19958    20157  Input           1  -
        4  |    chr1            10149    10348  Input           1  -
        3  |    chr1            10082    10281  Input           1  -
        2  |    chr1            10079    10278  Input           1  -
        0  |    chr1             9988    10187  Input           1  -
  PyRanges with 5 rows, 6 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.

  >>> gr2.loci['chr1', '+'] # get all rows with chromosome 'chr1' and strand '+'
    index  |    Chromosome      Start      End  Name        Score  Strand
    int64  |    category        int64    int64  object      int64  category
  -------  ---  ------------  -------  -------  --------  -------  ----------
        1  |    chr1            10073    10272  Input           1  +
        5  |    chr1            10280    10479  Input           1  +
        6  |    chr1            16056    16255  Input           1  +
        7  |    chr1            16064    16263  Input           1  +
        8  |    chr1            16109    16308  Input           1  +
  PyRanges with 5 rows, 6 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.

  >>> gr2.loci['chr1', 10000:11000] # get all rows on 'chr1' and overlapping 10000:11000
    index  |    Chromosome      Start      End  Name        Score  Strand
    int64  |    category        int64    int64  object      int64  category
  -------  ---  ------------  -------  -------  --------  -------  ----------
        1  |    chr1            10073    10272  Input           1  +
        5  |    chr1            10280    10479  Input           1  +
        4  |    chr1            10149    10348  Input           1  -
        3  |    chr1            10082    10281  Input           1  -
        2  |    chr1            10079    10278  Input           1  -
        0  |    chr1             9988    10187  Input           1  -
  PyRanges with 6 rows, 6 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

  >>> gr2.loci['chr1', '+', 10000:11000] # get all rows on 'chr1', strand '+', and overlapping 10000:11000
    index  |    Chromosome      Start      End  Name        Score  Strand
    int64  |    category        int64    int64  object      int64  category
  -------  ---  ------------  -------  -------  --------  -------  ----------
        1  |    chr1            10073    10272  Input           1  +
        5  |    chr1            10280    10479  Input           1  +
  PyRanges with 2 rows, 6 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.

  To use this kind of selection in combination with ``loc`` or ``iloc``, you can use the ``index`` attribute:

  >>> sindex=gr2.loci['chr1', '+', 10000:11000].index
  >>> gr2.loc[sindex, "Score"]=100
  >>> gr2
  index    |    Chromosome    Start    End      Name      Score    Strand
  int64    |    category      int64    int64    object    int64    category
  -------  ---  ------------  -------  -------  --------  -------  ----------
  1        |    chr1          10073    10272    Input     100      +
  5        |    chr1          10280    10479    Input     100      +
  6        |    chr1          16056    16255    Input     1        +
  7        |    chr1          16064    16263    Input     1        +
  ...      |    ...           ...      ...      ...       ...      ...
  4        |    chr1          10149    10348    Input     1        -
  3        |    chr1          10082    10281    Input     1        -
  2        |    chr1          10079    10278    Input     1        -
  0        |    chr1          9988     10187    Input     1        -
  PyRanges with 10 rows, 6 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

  Note that the above syntax works only for PyRanges with numerical indices.
  You can always call ``reset_index`` to make the index numerical.



From here!


Selecting columns
-----------------

As previously seen, single PyRanges column (which are pandas Series) can be extracted through the dot notation:


  >>> gr = pr.data.chipseq()
  >>> gr.Chromosome
  18      chr1
  70      chr1
  129     chr1
  170     chr1
  196     chr1
  	  ...
  3023    chrY
  3131    chrY
  3816    chrY
  3897    chrY
  9570    chrY
  Name: Chromosome, Length: 10000, dtype: category
  Categories (24, object): ['chr1', 'chr10', 'chr11', 'chr12', ..., 'chr8', 'chr9', 'chrX', 'chrY']

The same syntax can be used for the core PyRanges columns (Chromosome, Strand, Start, End) or for metadata columns:

  >>> gr.Name
  18      U0
  70      U0
  129     U0
  170     U0
  196     U0
  	  ..
  3023    U0
  3131    U0
  3816    U0
  3897    U0
  9570    U0
  Name: Name, Length: 10000, dtype: object

This syntax is analogous to pandas Dataframes. Note that, however, the bracket column selection in pandas does not work in the same way in PyRanges:

  >>> df=gr.df
  >>> df['End']
  0       212609559
  1       169887554
  2       216711036
  3       144227104
  4       148177850
            ...
  9995      7046834
  9996     15224260
  9997     13517917
  9998      8010976
  9999      7405401
  Name: End, Length: 10000, dtype: int64

  >>> gr['End']
  Empty PyRanges

Because the last expression is evaluated as a genomic region, i.e. a form of row selection: it is searching for intervals on a Chromosome named "End", and finds none. Indeed, this fetches intervals on the chrY:

  >>> gr['chrY']
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   | Start     | End       | Name       | Score     | Strand       |
  | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chrY         | 12930373  | 12930398  | U0         | 0         | +            |
  | chrY         | 15548022  | 15548047  | U0         | 0         | +            |
  | chrY         | 7194340   | 7194365   | U0         | 0         | +            |
  | chrY         | 21559181  | 21559206  | U0         | 0         | +            |
  | ...          | ...       | ...       | ...        | ...       | ...          |
  | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
  | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
  | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
  | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 23 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

You can provide a list of column names in the bracket notation to select those columns. Pyranges will still return a PyRanges object, therefore retaining the core columns regardless of whether they were selected or not:

  >>> gr[ ['Name'] ]
  +--------------+-----------+-----------+------------+--------------+
  | Chromosome   | Start     | End       | Name       | Strand       |
  | (category)   | (int64)   | (int64)   | (object)   | (category)   |
  |--------------+-----------+-----------+------------+--------------|
  | chr1         | 212609534 | 212609559 | U0         | +            |
  | chr1         | 169887529 | 169887554 | U0         | +            |
  | chr1         | 216711011 | 216711036 | U0         | +            |
  | chr1         | 144227079 | 144227104 | U0         | +            |
  | ...          | ...       | ...       | ...        | ...          |
  | chrY         | 15224235  | 15224260  | U0         | -            |
  | chrY         | 13517892  | 13517917  | U0         | -            |
  | chrY         | 8010951   | 8010976   | U0         | -            |
  | chrY         | 7405376   | 7405401   | U0         | -            |
  +--------------+-----------+-----------+------------+--------------+
  Stranded PyRanges object has 10,000 rows and 5 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

This is convenient to reduce genome annotation that consists of many columns:

  >>> ensembl_path = pr.get_example_path("ensembl.gtf")
  >>> ge = pr.read_gtf(ensembl_path)
  >>> ge
  +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+-----------------+----------------+-------------+----------------+-------+
  | Chromosome   | Source     | Feature      | Start     | End       | Score      | Strand       | Frame      | gene_id         | gene_version   | gene_name   | gene_source    | +14   |
  | (category)   | (object)   | (category)   | (int64)   | (int64)   | (object)   | (category)   | (object)   | (object)        | (object)       | (object)    | (object)       | ...   |
  |--------------+------------+--------------+-----------+-----------+------------+--------------+------------+-----------------+----------------+-------------+----------------+-------|
  | 1            | havana     | gene         | 11868     | 14409     | .          | +            | .          | ENSG00000223972 | 5              | DDX11L1     | havana         | ...   |
  | 1            | havana     | transcript   | 11868     | 14409     | .          | +            | .          | ENSG00000223972 | 5              | DDX11L1     | havana         | ...   |
  | 1            | havana     | exon         | 11868     | 12227     | .          | +            | .          | ENSG00000223972 | 5              | DDX11L1     | havana         | ...   |
  | 1            | havana     | exon         | 12612     | 12721     | .          | +            | .          | ENSG00000223972 | 5              | DDX11L1     | havana         | ...   |
  | ...          | ...        | ...          | ...       | ...       | ...        | ...          | ...        | ...             | ...            | ...         | ...            | ...   |
  | 1            | ensembl    | transcript   | 120724    | 133723    | .          | -            | .          | ENSG00000238009 | 6              | AL627309.1  | ensembl_havana | ...   |
  | 1            | ensembl    | exon         | 133373    | 133723    | .          | -            | .          | ENSG00000238009 | 6              | AL627309.1  | ensembl_havana | ...   |
  | 1            | ensembl    | exon         | 129054    | 129223    | .          | -            | .          | ENSG00000238009 | 6              | AL627309.1  | ensembl_havana | ...   |
  | 1            | ensembl    | exon         | 120873    | 120932    | .          | -            | .          | ENSG00000238009 | 6              | AL627309.1  | ensembl_havana | ...   |
  +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+-----------------+----------------+-------------+----------------+-------+
  Stranded PyRanges object has 95 rows and 26 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.
  14 hidden columns: gene_biotype, transcript_id, transcript_version, transcript_name, transcript_source, transcript_biotype, tag, transcript_support_level, exon_number, exon_id, exon_version, ... (+ 3 more.)


  >>> ge[ ['gene_id', 'gene_name'] ]
  +--------------+-----------+-----------+--------------+-----------------+-------------+
  | Chromosome   | Start     | End       | Strand       | gene_id         | gene_name   |
  | (category)   | (int64)   | (int64)   | (category)   | (object)        | (object)    |
  |--------------+-----------+-----------+--------------+-----------------+-------------|
  | 1            | 11868     | 14409     | +            | ENSG00000223972 | DDX11L1     |
  | 1            | 11868     | 14409     | +            | ENSG00000223972 | DDX11L1     |
  | 1            | 11868     | 12227     | +            | ENSG00000223972 | DDX11L1     |
  | 1            | 12612     | 12721     | +            | ENSG00000223972 | DDX11L1     |
  | ...          | ...       | ...       | ...          | ...             | ...         |
  | 1            | 120724    | 133723    | -            | ENSG00000238009 | AL627309.1  |
  | 1            | 133373    | 133723    | -            | ENSG00000238009 | AL627309.1  |
  | 1            | 129054    | 129223    | -            | ENSG00000238009 | AL627309.1  |
  | 1            | 120873    | 120932    | -            | ENSG00000238009 | AL627309.1  |
  +--------------+-----------+-----------+--------------+-----------------+-------------+
  Stranded PyRanges object has 95 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

The **drop method** is an alternative way of column selection wherein we specify what we want to remove, rather than what to keep:


  >>> gr.print()
  +--------------+-----------+-----------+------------+-----------+--------------+
  | Chromosome   | Start     | End       | Name       | Score     | Strand       |
  | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
  |--------------+-----------+-----------+------------+-----------+--------------|
  | chr1         | 212609534 | 212609559 | U0         | 0         | +            |
  | chr1         | 169887529 | 169887554 | U0         | 0         | +            |
  | chr1         | 216711011 | 216711036 | U0         | 0         | +            |
  | chr1         | 144227079 | 144227104 | U0         | 0         | +            |
  | ...          | ...       | ...       | ...        | ...       | ...          |
  | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
  | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
  | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
  | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
  +--------------+-----------+-----------+------------+-----------+--------------+
  Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  >>> gr.drop(['Name']).print()
  +--------------+-----------+-----------+-----------+--------------+
  | Chromosome   | Start     | End       | Score     | Strand       |
  | (category)   | (int64)   | (int64)   | (int64)   | (category)   |
  |--------------+-----------+-----------+-----------+--------------|
  | chr1         | 212609534 | 212609559 | 0         | +            |
  | chr1         | 169887529 | 169887554 | 0         | +            |
  | chr1         | 216711011 | 216711036 | 0         | +            |
  | chr1         | 144227079 | 144227104 | 0         | +            |
  | ...          | ...       | ...       | ...       | ...          |
  | chrY         | 15224235  | 15224260  | 0         | -            |
  | chrY         | 13517892  | 13517917  | 0         | -            |
  | chrY         | 8010951   | 8010976   | 0         | -            |
  | chrY         | 7405376   | 7405401   | 0         | -            |
  +--------------+-----------+-----------+-----------+--------------+
  Stranded PyRanges object has 10,000 rows and 5 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

Without arguments, drop will get rid of all non-core columns:

  >>> gr.drop()
  +--------------+-----------+-----------+--------------+
  | Chromosome   | Start     | End       | Strand       |
  | (category)   | (int64)   | (int64)   | (category)   |
  |--------------+-----------+-----------+--------------|
  | chr1         | 212609534 | 212609559 | +            |
  | chr1         | 169887529 | 169887554 | +            |
  | chr1         | 216711011 | 216711036 | +            |
  | chr1         | 144227079 | 144227104 | +            |
  | ...          | ...       | ...       | ...          |
  | chrY         | 15224235  | 15224260  | -            |
  | chrY         | 13517892  | 13517917  | -            |
  | chrY         | 8010951   | 8010976   | -            |
  | chrY         | 7405376   | 7405401   | -            |
  +--------------+-----------+-----------+--------------+
  Stranded PyRanges object has 10,000 rows and 4 columns from 24 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.


If you want to obtain a DataFrame with certain columns rather than a PyRanges object, get a DataFrame copy through the df property, then perform pandas-style column selection. Obviously, in this case core columns are returned only if explicitly selected:

  >>> gr.df [ ['Name', 'Start'] ]
       Name      Start
  0      U0  212609534
  1      U0  169887529
  2      U0  216711011
  3      U0  144227079
  4      U0  148177825
  ...   ...        ...
  9995   U0    7046809
  9996   U0   15224235
  9997   U0   13517892
  9998   U0    8010951
  9999   U0    7405376
  <BLANKLINE>
  [10000 rows x 2 columns]



Obtaining sequences
-------------------


A common operation is to fetch the sequences corresponding to the intervals represented in the PyRanges object. Function ``get_sequence`` takes as input a PyRanges object and the path to a fasta file, and returns a Series containing sequences, in the same order as the intervals. It requires package pyfaidx (install with pip install pyfaidx).

In the tutorial, we saw its usage with a real genome. Let's make a toy example here:

  >>> with open('minigenome.fa', 'w') as fw:
  ...     fw.write('>chrZ\n')
  ...     fw.write('AAAGGGCCCTTTAAAGGGCCCTTTAAAGGGCCCTTT\n')

  >>> sg = pr.from_dict({"Chromosome": ["chrZ", "chrZ", "chrZ", "chrZ"],
  ... 	           "Start": [0, 5, 10, 10], "End": [3, 8, 20, 20],
  ... 	           "name":["a", "a", "b", "c"],
  ... 	           "Strand":["+", "+", "+", "-"] })

  >>> sg
  +--------------+-----------+-----------+------------+--------------+
  | Chromosome   |     Start |       End | name       | Strand       |
  | (category)   |   (int64) |   (int64) | (object)   | (category)   |
  |--------------+-----------+-----------+------------+--------------|
  | chrZ         |         0 |         3 | a          | +            |
  | chrZ         |         5 |         8 | a          | +            |
  | chrZ         |        10 |        20 | b          | +            |
  | chrZ         |        10 |        20 | c          | -            |
  +--------------+-----------+-----------+------------+--------------+
  Stranded PyRanges object has 4 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

Note the genome sequence in the code above. Let's run ``get_sequences`` to obtain the portions corresponding to our intervals:


  >>> pr.get_sequence(sg, 'minigenome.fa')
  0           AAA
  1           GCC
  2    TTAAAGGGCC
  3    GGCCCTTTAA
  dtype: object

Note that the last two intervals have identical coordinates but are on opposite strands. Function ``get_sequence`` returns the reverse complement for intervals on the negative strand.

Since the returned Series has the same length as the PyRanges object, we can assign it to a new column:


  >>> sg.Sequence = pr.get_sequence(sg, 'minigenome.fa')
  >>> sg
  +--------------+-----------+-----------+------------+--------------+------------+
  | Chromosome   |     Start |       End | name       | Strand       | Sequence   |
  | (category)   |   (int64) |   (int64) | (object)   | (category)   | (object)   |
  |--------------+-----------+-----------+------------+--------------+------------|
  | chrZ         |         0 |         3 | a          | +            | AAA        |
  | chrZ         |         5 |         8 | a          | +            | GCC        |
  | chrZ         |        10 |        20 | b          | +            | TTAAAGGGCC |
  | chrZ         |        10 |        20 | c          | -            | GGCCCTTTAA |
  +--------------+-----------+-----------+------------+--------------+------------+
  Stranded PyRanges object has 4 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

This allows us to filter by sequence, using pandas string methods. For example, let's get those that start with G:



  >>> sg[sg.Sequence.str.startswith('G')]
  +--------------+-----------+-----------+------------+--------------+------------+
  | Chromosome   |     Start |       End | name       | Strand       | Sequence   |
  | (category)   |   (int64) |   (int64) | (object)   | (category)   | (object)   |
  |--------------+-----------+-----------+------------+--------------+------------|
  | chrZ         |         5 |         8 | a          | +            | GCC        |
  | chrZ         |        10 |        20 | c          | -            | GGCCCTTTAA |
  +--------------+-----------+-----------+------------+--------------+------------+
  Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

Let's get those which contain a CC and AA dinucleotides separated by 1-3 nucleotides:



  >>> sg[sg.Sequence.str.contains(r'CC.{1,3}AA', regex=True)]
  +--------------+-----------+-----------+------------+--------------+------------+
  | Chromosome   |     Start |       End | name       | Strand       | Sequence   |
  | (category)   |   (int64) |   (int64) | (object)   | (category)   | (object)   |
  |--------------+-----------+-----------+------------+--------------+------------|
  | chrZ         |        10 |        20 | c          | -            | GGCCCTTTAA |
  +--------------+-----------+-----------+------------+--------------+------------+
  Stranded PyRanges object has 1 rows and 6 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.



Function ``get_sequence`` will treat each interval independently. Often, you want to get the sequence of an mRNA, i.e. concatenating exons. Function get_transcript_sequence serves this purpose, and employs argument group_by to group the exons into mRNAs:


  >>> pr.get_transcript_sequence(sg, group_by='name', path='minigenome.fa')
    name    Sequence
  0    a      AAAGCC
  1    b  TTAAAGGGCC
  2    c  GGCCCTTTAA

Note that this returns a pandas DataFrame with a row per exon group: its shape is different from the original PyRanges.



Operating with data
~~~~~~~~~~~~~~~~~~~


In this section, we give an overview of methods to modify the data in PyRanges.
Changing row order
Methods sort allows to sort intervals, i.e. altering the order of rows in the PyRanges object. When run without arguments, orders interval by increasing Start. Commonly, genomic annotation files are sorted in this way.


  >>> sg = pr.from_dict({"Chromosome": ["chrA", "chrA", "chrB", "chrB", "chrB"],
  ... 	           "Start": [55, 20, 65, 35, 75],
  ... 	           "End": [88, 30, 75, 45, 85],
  ... 	           "name":["a", "a", "b", "c", "c"],
  ... 	           "Strand":["+", "+", "+", "-", "-"] })
  >>> sg
  +--------------+-----------+-----------+------------+--------------+
  | Chromosome   |     Start |       End | name       | Strand       |
  | (category)   |   (int64) |   (int64) | (object)   | (category)   |
  |--------------+-----------+-----------+------------+--------------|
  | chrA         |        55 |        88 | a          | +            |
  | chrA         |        20 |        30 | a          | +            |
  | chrB         |        65 |        75 | b          | +            |
  | chrB         |        35 |        45 | c          | -            |
  | chrB         |        75 |        85 | c          | -            |
  +--------------+-----------+-----------+------------+--------------+
  Stranded PyRanges object has 5 rows and 5 columns from 2 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  >>> sg.sort()
  +--------------+-----------+-----------+------------+--------------+
  | Chromosome   |     Start |       End | name       | Strand       |
  | (category)   |   (int64) |   (int64) | (object)   | (category)   |
  |--------------+-----------+-----------+------------+--------------|
  | chrA         |        20 |        30 | a          | +            |
  | chrA         |        55 |        88 | a          | +            |
  | chrB         |        65 |        75 | b          | +            |
  | chrB         |        35 |        45 | c          | -            |
  | chrB         |        75 |        85 | c          | -            |
  +--------------+-----------+-----------+------------+--------------+
  Stranded PyRanges object has 5 rows and 5 columns from 2 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.


Remember that **sorting is performed separately for each internal table**: intervals on different chromosome/strands won't ever cross each other. To have all intervals sorted, work with a DataFrame object instead.

For intervals on the negative strand, it may be convenient to sort in the opposite order, since for them the leftmost exon is actually the last one in the mRNA. Instead of having to split the PyRanges object for this task, you may run sort with special argument "5", which will sort intervals in 5' to 3' order:


  >>> sg.sort('5')
  +--------------+-----------+-----------+------------+--------------+
  | Chromosome   |     Start |       End | name       | Strand       |
  | (category)   |   (int64) |   (int64) | (object)   | (category)   |
  |--------------+-----------+-----------+------------+--------------|
  | chrA         |        20 |        30 | a          | +            |
  | chrA         |        55 |        88 | a          | +            |
  | chrB         |        65 |        75 | b          | +            |
  | chrB         |        75 |        85 | c          | -            |
  | chrB         |        35 |        45 | c          | -            |
  +--------------+-----------+-----------+------------+--------------+
  Stranded PyRanges object has 5 rows and 5 columns from 2 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

Sorting may also take any column name, or a list of colum names, to sort rows by their value:

  >>> ag = pr.from_dict({"Chromosome": "chrX",
  ... 	           "Start": [55, 65, 20, 35, 75],
  ... 	           "End": [88, 75, 30, 45, 85],
  ... 	           "Strand":["+", "+", "+", "+", "+"],
  ... 	           "col1":[1, 4, 4, 2, 2],
  ... 	            })
  >>> ag
  +--------------+-----------+-----------+--------------+-----------+
  | Chromosome   |     Start |       End | Strand       |      col1 |
  | (category)   |   (int64) |   (int64) | (category)   |   (int64) |
  |--------------+-----------+-----------+--------------+-----------|
  | chrX         |        55 |        88 | +            |         1 |
  | chrX         |        65 |        75 | +            |         4 |
  | chrX         |        20 |        30 | +            |         4 |
  | chrX         |        35 |        45 | +            |         2 |
  | chrX         |        75 |        85 | +            |         2 |
  +--------------+-----------+-----------+--------------+-----------+
  Stranded PyRanges object has 5 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  >>> ag.sort('col1')
  +--------------+-----------+-----------+--------------+-----------+
  | Chromosome   |     Start |       End | Strand       |      col1 |
  | (category)   |   (int64) |   (int64) | (category)   |   (int64) |
  |--------------+-----------+-----------+--------------+-----------|
  | chrX         |        55 |        88 | +            |         1 |
  | chrX         |        35 |        45 | +            |         2 |
  | chrX         |        75 |        85 | +            |         2 |
  | chrX         |        65 |        75 | +            |         4 |
  | chrX         |        20 |        30 | +            |         4 |
  +--------------+-----------+-----------+--------------+-----------+
  Stranded PyRanges object has 5 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.

  >>> ag.sort(['col1', 'End'])
  +--------------+-----------+-----------+--------------+-----------+
  | Chromosome   |     Start |       End | Strand       |      col1 |
  | (category)   |   (int64) |   (int64) | (category)   |   (int64) |
  |--------------+-----------+-----------+--------------+-----------|
  | chrX         |        55 |        88 | +            |         1 |
  | chrX         |        35 |        45 | +            |         2 |
  | chrX         |        75 |        85 | +            |         2 |
  | chrX         |        20 |        30 | +            |         4 |
  | chrX         |        65 |        75 | +            |         4 |
  +--------------+-----------+-----------+--------------+-----------+
  Stranded PyRanges object has 5 rows and 5 columns from 1 chromosomes.
  For printing, the PyRanges was sorted on Chromosome and Strand.


..
    [add note: index are not allowed. Stil, you can use sort to get rows in a certain order]
    Operations on coordinates
    [change columns as series: p.Start+=1000 ...]
    [... however there are more convenient methods: subsequence, spliced_sequence, extend]
    [after extend, show genome_bounds]

    Operations on metadata columns:
    [insert new columns: 1. p.Col1=... or 2. assign method. 3. Assign with multiple ones at once]

    Operations on multiple pyranges
    [concatenation: use pandas and turn to pyranges]

    A common operation on (multiple) pyranges regard overlaps. These are shown in the next page


    Overlapping and matching PyRanges
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    [present different methods for different aims that all have to do with overlap: merge, cluster, subtract, join, count_overlaps ... . Start with a table summarizing differences: input, output].
    [add note: pandas merge: different!]

    Summary statistics
    ~~~~~~~~~~~~~~~~~~

    [Create count-matrix from multiple PyRanges]
    [all stats methods presented briefly]

    Computing with PyRanges
    ~~~~~~~~~~~~~~~~~~~~~~~

    [ready made methods should cover most things]
    [possibility to chain things to save memory]
    [outline strategies for custom methods: apply and similar methods]
    [Also cite the simple but not optimal: convert to dataframes / or iterate through groups of same-chrom dataframes]
    [multiple cores]

    Working at the transcript level
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    [spliced_subsequence, subsequence, get_transcript_sequence,
    extend (to be developed with group_by),
    boundaries ,
    cumsum groupby as example

    ]


    Fetching external gene tracks
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    [if pyranges_db is a thing, describe its uses here]


    RLEs: run length encodings
    ~~~~~~~~~~~~~~~~~~~~~~~~~~

    [outline as advanced usage. Put everything related to RLEs in a single chapter; keep as last even if you add further chapters]


