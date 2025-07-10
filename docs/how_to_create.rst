Loading/Creating PyRanges
~~~~~~~~~~~~~~~~~~~~~~~~~

.. contents::
   :local:
   :depth: 2


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

You can instantiate a PyRanges object using a dictionary:

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


As in the creation of a Dataframe, each list in the dictionary must have the same length, i.e. the number of rows
in the PyRanges. Also, it may be any iterable, including pd.Series or np.array.
Alternatively, if a string or scalar is provided, it is broadcasted to the length of the other columns:

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


To read bam files, the optional bamread-library must be installed, with ::

    pip install bamread

Let's read a bam file:

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

``read_bam`` takes various arguments, such as ``sparse``.
With ``sparse=True`` (default), only the columns ``['Chromosome', 'Start', 'End', 'Strand', 'Flag']``
are fetched. Setting ``sparse=False`` additionally gives you the columns
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

Concatenating PyRanges
----------------------

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

Data for testing
----------------

For testing purposes, pyranges provides some data in :ref:`pr.example_data <pyranges_example_data>`.
See an overview with:

  >>> pr.example_data
  Available example data:
  -----------------------
  example_data.chipseq            : Example ChIP-seq data.
  example_data.chipseq_background : Example ChIP-seq data.
  example_data.chromsizes         : Example chromsizes data (hg19).
  example_data.ensembl_gtf        : Example gtf file from Ensembl.
  example_data.interpro_hits      : Example of InterPro protein hits.
  example_data.rfam_hits          : Example of RNA motifs (Rfam) as 1-based dataframe.
  example_data.f1                 : Example bed file.
  example_data.f2                 : Example bed file.
  example_data.aorta              : Example ChIP-seq data.
  example_data.aorta2             : Example ChIP-seq data.
  example_data.ncbi_gff           : Example NCBI GFF data.
  example_data.ncbi_fasta         : Example NCBI fasta.
  example_data.files              : Return a dict of basenames to file paths of available data.

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
        0  |    chr12         108700348  108700398  +
        1  |    chr1          230144267  230144317  -
        2  |    chr3           54767920   54767970  +
        3  |    chr3          162329749  162329799  +
        4  |    chr3          176218669  176218719  +
  PyRanges with 5 rows, 4 columns, and 1 index columns.
  Contains 3 chromosomes and 2 strands.
