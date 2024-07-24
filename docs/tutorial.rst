Tutorial
========

:class:`PyRanges <pyranges.PyRanges>` objects represent sequence intervals ("ranges") such as genomic regions.
Examples include gene annotations, mapped reads, protein domains, and results from
Blast or analogous software. Pyranges provides convenient methods to load data in
GFF, GTF, BED, BAM format. In this tutorial, we refer to PyRanges objects for their
typical use case, i.e. representing genomic intervals, but the same concepts and methods
may be applied to RNA or protein sequences.

Every interval in a PyRanges object is defined, at minimum, by its chromosome and start
and end coordinates. Optionally, the strand (+ or -) can also be present; if so, the
PyRanges object is "Stranded". The ranges can additionally have an arbitrary number
of meta-data fields, i.e. columns associated with them.

Note that PyRanges follows the standard python convention:

* coordinates are **0-based**
* in each interval, the **start is included** and the **end is excluded**.

Some genomic coordinate formats (GFF, GTF) use a 1-based, start-and-end-included format.
Pyranges takes care of converting between these conventions when loading and writing files in these formats.

PyRanges objects are a subclass of  `Pandas <https://pandas.pydata.org/>`_ DataFrame, extending its functionalities
to genomic operations. This means that the vast pandas/numpy ecosystem for high-performance scientific computing is
available to manipulate the data directly in PyRanges objects. While not strictly necessary, having
familiarity with pandas greatly facilitates the use of PyRanges.

.. contents:: Contents of Tutorial
   :depth: 3


Getting started
~~~~~~~~~~~~~~~

We recommend using `ipython <https://ipython.readthedocs.io/>`_ or `Jupyter <https://jupyter.org/>`_ for this tutorial.
Besides pyranges and pandas, we will use optional modules (e.g. **pyfaidx**).
If you haven't already, install the optional add-ons with:

.. code-block:: shell

      pip install pyranges1[add-ons]


Loading and accessing pyranges objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's import libraries, and load in memory an example annotation in GFF3 format, consisting of a portion of the genome
annotation of the worm *Dimorphilus gyrociliatus*.

  >>> import pyranges as pr
  >>> ann = pr.example_data.ncbi_gff
  >>> ann
  index    |    Chromosome         Source    Feature     Start    End      Score     Strand      Frame     ...
  int64    |    category           object    category    int64    int64    object    category    object    ...
  -------  ---  -----------------  --------  ----------  -------  -------  --------  ----------  --------  -----
  0        |    CAJFCJ010000053.1  EMBL      region      0        109277   .         +           .         ...
  1        |    CAJFCJ010000053.1  EMBL      gene        4882     5264     .         -           .         ...
  2        |    CAJFCJ010000053.1  EMBL      mRNA        4882     5264     .         -           .         ...
  3        |    CAJFCJ010000053.1  EMBL      exon        4882     5264     .         -           .         ...
  ...      |    ...                ...       ...         ...      ...      ...       ...         ...       ...
  146      |    CAJFCJ010000025.1  EMBL      CDS         2753     2851     .         -           0         ...
  147      |    CAJFCJ010000025.1  EMBL      CDS         2593     2693     .         -           1         ...
  148      |    CAJFCJ010000025.1  EMBL      CDS         2354     2537     .         -           0         ...
  149      |    CAJFCJ010000025.1  EMBL      CDS         2174     2294     .         -           0         ...
  PyRanges with 150 rows, 21 columns, and 1 index columns. (13 columns not shown: "ID", "Dbxref", "gbkey", ...).
  Contains 6 chromosomes and 2 strands.

Above, we leveraged the builtin example data. In real use cases, you would load data from a file, using methods like
:func:`read_gtf <pyranges.read_gtf>`, :func:`read_bed <pyranges.read_bed>` and others (see all readers in
:doc:`pyranges_module`).

The ``ann`` object has lots of columns, most of which are not displayed because of space. Here's the full list:

  >>> ann.columns
  Index(['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand',
         'Frame', 'ID', 'Dbxref', 'gbkey', 'mol_type', 'note', 'Name',
         'gene_biotype', 'locus_tag', 'Parent', 'Note', 'standard_name',
         'product', 'protein_id'],
        dtype='object')


Let's select only certain columns. We can use the method
:func:`get_with_loc_columns <pyranges.PyRanges.get_with_loc_columns>` to select columns by name, and
retain the "genomic location" columns **Chromosome, Start, End**, (and **Strand** if present):

  >>> ann = ann.get_with_loc_columns(['Feature', 'Parent', 'ID'])
  >>> ann
  index    |    Chromosome         Start    End      Strand      Feature     Parent                 ...
  int64    |    category           int64    int64    category    category    object                 ...
  -------  ---  -----------------  -------  -------  ----------  ----------  ---------------------  -----
  0        |    CAJFCJ010000053.1  0        109277   +           region      nan                    ...
  1        |    CAJFCJ010000053.1  4882     5264     -           gene        nan                    ...
  2        |    CAJFCJ010000053.1  4882     5264     -           mRNA        gene-DGYR_LOCUS13733   ...
  3        |    CAJFCJ010000053.1  4882     5264     -           exon        rna-DGYR_LOCUS13733    ...
  ...      |    ...                ...      ...      ...         ...         ...                    ...
  146      |    CAJFCJ010000025.1  2753     2851     -           CDS         rna-DGYR_LOCUS12552-2  ...
  147      |    CAJFCJ010000025.1  2593     2693     -           CDS         rna-DGYR_LOCUS12552-2  ...
  148      |    CAJFCJ010000025.1  2354     2537     -           CDS         rna-DGYR_LOCUS12552-2  ...
  149      |    CAJFCJ010000025.1  2174     2294     -           CDS         rna-DGYR_LOCUS12552-2  ...
  PyRanges with 150 rows, 7 columns, and 1 index columns. (1 columns not shown: "ID").
  Contains 6 chromosomes and 2 strands.

The Chromosome column can take any value among the sequence names in the genome assembly.
In top-quality assemblies, it corresponds to actual chromosomes, and in other cases it is contigs or scaffolds;
for simplicity, here we refer to it as chromosomes. In a fasta file, the sequence name is the first word of a header
line (i.e. those starting with ">"). Let's peek the assembly fasta file available as example data:

  >>> genome_file = pr.example_data.files['ncbi.fasta']
  >>> with open(genome_file) as fh:
  ...   for _ in range(8):
  ...     print(fh.readline().strip())
  >CAJFCJ010000053.1 Dimorphilus gyrociliatus genome assembly, contig: scaffold053, whole genome shotgun sequence
  aaaaaaagaagtttttgacaaactttttctttttttcatcaagCTTTGTATAATGGACAA
  ACTAACgcaactttttcaattactGTTAACAAACTACCTGAAACAATTTACAATTCAAAA
  AGTACATTTTGTATTAGAAATTATTCCAAGAAAATTCAAGTAGATTTGAAATTCATGATT
  TAACTTGTGAAATTGTGTataaggaaaatatataaatattttcaaaactgTTACTTTGGA
  TACTAAAGAAATTCCattagaaataattgaaatatttgtatatacttcaccaaatgaaag
  aatgaatgaaataagtaaaaataaaatggagaaatttttttttttaattttttttctctt
  tcttcctttattCATAGctttatttgataatttcaaGAGTATAATTGAAGAGATCAGTGT


Genomic annotations often contain information for diverse entities, such as genes, mRNAs, exons, CDS, etc.
In GFF files, the entity type is encoded in the Feature column. In pyranges, you use the dot notation to
fetch an individual column, which is technically a pandas Series:

  >>> ann.Feature # or ann['Feature']
  0      region
  1        gene
  2        mRNA
  3        exon
  4         CDS
          ...
  145       CDS
  146       CDS
  147       CDS
  148       CDS
  149       CDS
  Name: Feature, Length: 150, dtype: category
  Categories (5, object): ['CDS', 'exon', 'gene', 'mRNA', 'region']


The syntax ``ann[column_name]`` is also available, and must be used when creating or updating a column.
Let's create a new column with the midpoint of each interval:

  >>> ann['midpoint'] = (ann.Start + ann.End) // 2
  >>> ann.get_with_loc_columns(['midpoint'])
  index    |    Chromosome         Start    End      Strand      midpoint
  int64    |    category           int64    int64    category    int64
  -------  ---  -----------------  -------  -------  ----------  ----------
  0        |    CAJFCJ010000053.1  0        109277   +           54638
  1        |    CAJFCJ010000053.1  4882     5264     -           5073
  2        |    CAJFCJ010000053.1  4882     5264     -           5073
  3        |    CAJFCJ010000053.1  4882     5264     -           5073
  ...      |    ...                ...      ...      ...         ...
  146      |    CAJFCJ010000025.1  2753     2851     -           2802
  147      |    CAJFCJ010000025.1  2593     2693     -           2643
  148      |    CAJFCJ010000025.1  2354     2537     -           2445
  149      |    CAJFCJ010000025.1  2174     2294     -           2234
  PyRanges with 150 rows, 5 columns, and 1 index columns.
  Contains 6 chromosomes and 2 strands.

Let's focus on a row subset of the annotation: CDS intervals, corresponding to coding sequences.
We filter rows and create a new PyRanges object called ``cds``:

  >>> selector = (ann.Feature == 'CDS')
  >>> cds = ann [selector]

The object ``selector`` is a Series of boolean values, so it can be used to index PyRanges.

Now, let's further reduce the width of the cds object.
We showcase an alternative method for column selection: ``drop`` lets us choose which columns to discard.

  >>> cds = cds.drop( ['Feature', 'Parent', 'midpoint'], axis=1 )
  >>> cds
  index    |    Chromosome         Start    End      Strand      ID
  int64    |    category           int64    int64    category    object
  -------  ---  -----------------  -------  -------  ----------  ----------------
  4        |    CAJFCJ010000053.1  4882     5263     -           cds-CAD5126491.1
  11       |    CAJFCJ010000053.1  10732    10958    +           cds-CAD5126492.1
  12       |    CAJFCJ010000053.1  11028    11169    +           cds-CAD5126492.1
  13       |    CAJFCJ010000053.1  11227    11400    +           cds-CAD5126492.1
  ...      |    ...                ...      ...      ...         ...
  146      |    CAJFCJ010000025.1  2753     2851     -           cds-CAD5125114.1
  147      |    CAJFCJ010000025.1  2593     2693     -           cds-CAD5125114.1
  148      |    CAJFCJ010000025.1  2354     2537     -           cds-CAD5125114.1
  149      |    CAJFCJ010000025.1  2174     2294     -           cds-CAD5125114.1
  PyRanges with 56 rows, 5 columns, and 1 index columns.
  Contains 3 chromosomes and 2 strands.


``drop`` is actually a method of pandas dataframe, inherited by PyRanges.
Whenever a pandas methods is applied to a PyRanges object, if the returned object has the genomic location columns,
then it is returned as a PyRanges object. Otherwise, a dataframe is returned.

We already seen a boolean selector to filter rows. The ``loc`` and ``iloc`` pandas operators are also available.
Besides, pyranges offers the :func:`loci <pyranges.PyRanges.loci>` operator for selecting intervals in a
genomic region of interest. It accepts various syntaxes.
The code below will show intervals overlapping with the specified position range in the requested chromosome:

  >>> reg = cds.loci['CAJFCJ010000097.1', '+', 50000:55000]
  >>> reg
  index    |    Chromosome         Start    End      Strand      ID
  int64    |    category           int64    int64    category    object
  -------  ---  -----------------  -------  -------  ----------  ----------------
  110      |    CAJFCJ010000097.1  51865    52382    +           cds-CAD5126878.1
  111      |    CAJFCJ010000097.1  52446    52826    +           cds-CAD5126878.1
  112      |    CAJFCJ010000097.1  52903    53027    +           cds-CAD5126878.1
  113      |    CAJFCJ010000097.1  53339    53404    +           cds-CAD5126878.1
  ...      |    ...                ...      ...      ...         ...
  121      |    CAJFCJ010000097.1  52261    52382    +           cds-CAD5126877.1
  122      |    CAJFCJ010000097.1  52446    52826    +           cds-CAD5126877.1
  123      |    CAJFCJ010000097.1  52903    53027    +           cds-CAD5126877.1
  124      |    CAJFCJ010000097.1  53339    53404    +           cds-CAD5126877.1
  PyRanges with 9 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.

We cannot see all rows because of space. We can set how many rows are displayed using
:func:`pyranges.options.set_option`:

  >>> pr.options.set_option('max_rows_to_show', 10)
  >>> reg
    index  |    Chromosome           Start      End  Strand      ID
    int64  |    category             int64    int64  category    object
  -------  ---  -----------------  -------  -------  ----------  ----------------
      110  |    CAJFCJ010000097.1    51865    52382  +           cds-CAD5126878.1
      111  |    CAJFCJ010000097.1    52446    52826  +           cds-CAD5126878.1
      112  |    CAJFCJ010000097.1    52903    53027  +           cds-CAD5126878.1
      113  |    CAJFCJ010000097.1    53339    53404  +           cds-CAD5126878.1
      120  |    CAJFCJ010000097.1    51865    52201  +           cds-CAD5126877.1
      121  |    CAJFCJ010000097.1    52261    52382  +           cds-CAD5126877.1
      122  |    CAJFCJ010000097.1    52446    52826  +           cds-CAD5126877.1
      123  |    CAJFCJ010000097.1    52903    53027  +           cds-CAD5126877.1
      124  |    CAJFCJ010000097.1    53339    53404  +           cds-CAD5126877.1
  PyRanges with 9 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.

Let's go back to default display settings:

  >>> pr.options.reset_options()

Working with groups of exons
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multi-exonic genes are represented with multiple rows in PyRanges. In this tutorial, the ``ID`` column links the
intervals belonging to the same CDS: these rows have the same ID value.
While this concept applies to all annotations, files from different sources may use different column names
for this purpose (e.g. transcript_id). Note that here we focus on CDS regions. These may encompass multiple exons,
but they do not span the whole mRNA: the 5'UTRs and 3'UTRs are not included.
Various PyRanges methods are available to work with groups of intervals, accepting argument ``transcript_id``.

Next, we will examine the first and last codon of annotated CDSs.
We will obtain their genomic coordinate, then fetch their sequence.

Method :func:`spliced_subsequence <pyranges.PyRanges.spliced_subsequence>` allows to obtain a subregion of
groups of intervals. The code below derives the first codon of each CDS group; grouping is defined by their ID:

  >>> first=cds.spliced_subsequence(start=0, end=3, transcript_id='ID')
  >>> first
  index    |    Chromosome         Start    End      Strand      ID
  int64    |    category           int64    int64    category    object
  -------  ---  -----------------  -------  -------  ----------  ----------------
  4        |    CAJFCJ010000053.1  5260     5263     -           cds-CAD5126491.1
  11       |    CAJFCJ010000053.1  10732    10735    +           cds-CAD5126492.1
  18       |    CAJFCJ010000053.1  19649    19652    +           cds-CAD5126493.1
  25       |    CAJFCJ010000053.1  27136    27139    -           cds-CAD5126494.1
  ...      |    ...                ...      ...      ...         ...
  120      |    CAJFCJ010000097.1  51865    51868    +           cds-CAD5126877.1
  135      |    CAJFCJ010000025.1  2753     2755     -           cds-CAD5125115.1
  136      |    CAJFCJ010000025.1  2692     2693     -           cds-CAD5125115.1
  145      |    CAJFCJ010000025.1  3150     3153     -           cds-CAD5125114.1
  PyRanges with 18 rows, 5 columns, and 1 index columns.
  Contains 3 chromosomes and 2 strands.

Let's **fetch the sequence** for each of these intervals from our genome fasta file.
The function :func:`get_sequence <pyranges.PyRanges.get_sequence>` returns one sequence per interval, which we assign to a new column of our pyranges object:

  >>> first['Sequence'] = first.get_sequence(genome_file)  #genome_file defined above
  >>> first
  index    |    Chromosome         Start    End      Strand      ID                Sequence
  int64    |    category           int64    int64    category    object            object
  -------  ---  -----------------  -------  -------  ----------  ----------------  ----------
  4        |    CAJFCJ010000053.1  5260     5263     -           cds-CAD5126491.1  ATG
  11       |    CAJFCJ010000053.1  10732    10735    +           cds-CAD5126492.1  ATG
  18       |    CAJFCJ010000053.1  19649    19652    +           cds-CAD5126493.1  ATG
  25       |    CAJFCJ010000053.1  27136    27139    -           cds-CAD5126494.1  ATG
  ...      |    ...                ...      ...      ...         ...               ...
  120      |    CAJFCJ010000097.1  51865    51868    +           cds-CAD5126877.1  ATG
  135      |    CAJFCJ010000025.1  2753     2755     -           cds-CAD5125115.1  at
  136      |    CAJFCJ010000025.1  2692     2693     -           cds-CAD5125115.1  g
  145      |    CAJFCJ010000025.1  3150     3153     -           cds-CAD5125114.1  ATG
  PyRanges with 18 rows, 6 columns, and 1 index columns.
  Contains 3 chromosomes and 2 strands.



The ``Sequence`` column is a pandas Series containing strings. We see that the starting codon is ATG in most cases, as expected.
When we check the length of the sequences, we notice that some are not 3-letter long:

  >>> bool( (first.Sequence.str.len() == 3 ).all() )
  False

Let's look at those sequences, using a row selector as before:

  >>> first [ first.Sequence.str.len() != 3 ]
    index  |    Chromosome           Start      End  Strand      ID                Sequence
    int64  |    category             int64    int64  category    object            object
  -------  ---  -----------------  -------  -------  ----------  ----------------  ----------
      135  |    CAJFCJ010000025.1     2753     2755  -           cds-CAD5125115.1  at
      136  |    CAJFCJ010000025.1     2692     2693  -           cds-CAD5125115.1  g
  PyRanges with 2 rows, 6 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.
  

In some cases the starting codon is split between two exons. This is uncommon, but expected at least in a few genes
in a genome. How do we get the full codon sequence?

Instead of :func:`get_sequence <pyranges.PyRanges.get_sequence>`, let's use
:func:`get_transcript_sequence <pyranges.PyRanges.get_transcript_sequence>` ,
which returns the concatenated sequence of a group of intervals,
i.e. joining exons together. The sequence is given 5' to 3'.

  >>> seq_first = first.get_transcript_sequence(transcript_id='ID', path=genome_file)
  >>> seq_first
                    ID Sequence
  0   cds-CAD5125114.1      ATG
  1   cds-CAD5125115.1      atg
  2   cds-CAD5126491.1      ATG
  3   cds-CAD5126492.1      ATG
  4   cds-CAD5126493.1      ATG
  5   cds-CAD5126494.1      ATG
  6   cds-CAD5126495.1      ATG
  7   cds-CAD5126496.1      atg
  8   cds-CAD5126497.1      ATG
  9   cds-CAD5126498.1      atg
  10  cds-CAD5126499.1      atg
  11  cds-CAD5126873.1      ATG
  12  cds-CAD5126874.1      ATG
  13  cds-CAD5126875.1      ATG
  14  cds-CAD5126876.1      ATG
  15  cds-CAD5126877.1      ATG
  16  cds-CAD5126878.1      ATG


``seq_first`` is not a PyRanges object, but a pandas DataFrame. It has a column for the group (ID) and one for Sequence.
Here we confirm the sequence length is always 3:

  >>> bool( (seq_first.Sequence.str.len()==3).all() )
  True


Ok, so far we got the coordinates and sequences of the first codon of each CDS.

Now let's look at  stop codons.
First, we get the a pyranges object of the last codon of each CDS.
Conveniently, :func:`spliced_subsequence <pyranges.PyRanges.spliced_subsequence>` accepts negative arguments
to count from the 3', so we can obtain the last three nucleotides of CDSs with:

  >>> last = cds.spliced_subsequence(start=-3, transcript_id='ID')

By not providing an ``end`` argument, we requested intervals that reach the very end of each CDS group.
Let's get their sequence as before:

  >>> seq_last = last.get_transcript_sequence(transcript_id='ID', path=genome_file)
  >>> seq_last['Sequence'] = seq_last['Sequence'].str.upper()
  >>> seq_last
                    ID Sequence
  0   cds-CAD5125114.1      TGA
  1   cds-CAD5125115.1      TGA
  2   cds-CAD5126491.1      TAA
  3   cds-CAD5126492.1      TGA
  4   cds-CAD5126493.1      TAA
  5   cds-CAD5126494.1      TAG
  6   cds-CAD5126495.1      TAA
  7   cds-CAD5126496.1      TGA
  8   cds-CAD5126497.1      TAA
  9   cds-CAD5126498.1      TAA
  10  cds-CAD5126499.1      TAG
  11  cds-CAD5126873.1      TGA
  12  cds-CAD5126874.1      TAG
  13  cds-CAD5126875.1      TAA
  14  cds-CAD5126876.1      TGA
  15  cds-CAD5126877.1      TAA
  16  cds-CAD5126878.1      TAA


Let's use pandas ``value_counts`` to see the usage of stop codons:

  >>> seq_last['Sequence'].value_counts()
  Sequence
  TAA    8
  TGA    6
  TAG    3
  Name: count, dtype: int64

Say we want to focus on CDSs with a TAA stop codon. Let's gather the IDs of those CDSs:

  >>> taa_stop_ids = seq_last[ seq_last.Sequence == 'TAA' ].ID

We can now use this list to subset the ``cds`` object:

  >>> taa_stop_cds = cds[ cds.ID.isin(taa_stop_ids) ]


Writing coordinates and sequences to the disk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We obtained a custom genome annotation, consisting of CDS with a TAA stop codon.
We can now write this :class:`PyRanges <pyranges.PyRanges>`
object to a file, for example in GTF format:

  >>> taa_stop_cds.to_gtf('Dgyro.taa_CDS.gtf')


Let's get the sequence for these CDSs and write it to a tabular file using pandas method ``to_csv``:

  >>> taa_stop_cds_seqs = taa_stop_cds.get_transcript_sequence(transcript_id='ID', path=genome_file)
  >>> taa_stop_cds_seqs.to_csv('Dgyro_taa_CDS_seqs.tsv', sep='\t', index=False)

Note that ``taa_stop_cds_seqs`` is a pandas DataFrame. To write sequences in fasta format we use:

  >>> with open('Dgyro_taa_CDS_seqs.fa', 'w') as fw: # doctest: +SKIP
  ...   for xin, xid, xseq in taa_stop_cds_seqs.itertuples():
  ...     fw.write(f'>{xid}\n{xseq}\n')


Extending genomic intervals
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now we want to obtain (a practical approximation of) promoter sequences, here defined as the
300bp region before the start codon. Before we begin, let's peek into our object ``cds`` using
the pandas method ``head``:

  >>> cds.head()
    index  |    Chromosome           Start      End  Strand      ID
    int64  |    category             int64    int64  category    object
  -------  ---  -----------------  -------  -------  ----------  ----------------
        4  |    CAJFCJ010000053.1     4882     5263  -           cds-CAD5126491.1
       11  |    CAJFCJ010000053.1    10732    10958  +           cds-CAD5126492.1
       12  |    CAJFCJ010000053.1    11028    11169  +           cds-CAD5126492.1
       13  |    CAJFCJ010000053.1    11227    11400  +           cds-CAD5126492.1
       14  |    CAJFCJ010000053.1    11453    14183  +           cds-CAD5126492.1
  PyRanges with 5 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

First, we use the method  :func:`extend <pyranges.PyRanges.extend>`
to obtain intervals which include the CDS and the promoter defined as above.
We will group by the ID column, so that the extension is applied to each CDS group
(i.e. in this case only the 5' most
interval of each group).

  >>> g = cds.extend(ext_5=300, transcript_id='ID')
  >>> g.head()
    index  |    Chromosome           Start      End  Strand      ID
    int64  |    category             int64    int64  category    object
  -------  ---  -----------------  -------  -------  ----------  ----------------
        0  |    CAJFCJ010000025.1     2753     3055  -           cds-CAD5125115.1
        1  |    CAJFCJ010000025.1     2593     2693  -           cds-CAD5125115.1
        2  |    CAJFCJ010000025.1     2354     2537  -           cds-CAD5125115.1
        3  |    CAJFCJ010000025.1     2174     2294  -           cds-CAD5125115.1
        4  |    CAJFCJ010000025.1     3111     3453  -           cds-CAD5125114.1
  PyRanges with 5 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.


In the object we obtained, the promoter corresponds to the first 300 bp of every interval group.
We can use method :func:`spliced_subsequence <pyranges.PyRanges.spliced_subsequence>`  again to get it:

  >>> prom = g.spliced_subsequence(0, 300, transcript_id='ID')
  >>> prom.head()
    index  |    Chromosome           Start      End  Strand      ID
    int64  |    category             int64    int64  category    object
  -------  ---  -----------------  -------  -------  ----------  ----------------
        0  |    CAJFCJ010000025.1     2755     3055  -           cds-CAD5125115.1
        4  |    CAJFCJ010000025.1     3153     3453  -           cds-CAD5125114.1
        9  |    CAJFCJ010000053.1    10432    10732  +           cds-CAD5126492.1
       13  |    CAJFCJ010000053.1    19349    19649  +           cds-CAD5126493.1
       14  |    CAJFCJ010000053.1    38860    39160  +           cds-CAD5126495.1
  PyRanges with 5 rows, 5 columns, and 1 index columns.
  Contains 2 chromosomes and 2 strands.


Because we extended intervals, some may have gone out-of-bounds on the left or on the right side:
they may have a Start smaller than 0, or an End greater than the length of its chromosome, respectively.
The function :func:`genome_bounds <pyranges.PyRanges.genome_bounds>`
is designed to correct this.
We may use it to remove out-of-bounds intervals, or to retain only their in-bound portions.
We go for the second option, with ``clip=True``:

  >>> import pyfaidx
  >>> pyf=pyfaidx.Fasta(genome_file)
  >>> cor_prom = prom.genome_bounds(chromsizes=pyf, clip=True)
  >>> cor_prom.head()
    index  |    Chromosome           Start      End  Strand      ID
    int64  |    category             int64    int64  category    object
  -------  ---  -----------------  -------  -------  ----------  ----------------
        0  |    CAJFCJ010000025.1     2755     3055  -           cds-CAD5125115.1
        4  |    CAJFCJ010000025.1     3153     3418  -           cds-CAD5125114.1
        9  |    CAJFCJ010000053.1    10432    10732  +           cds-CAD5126492.1
       13  |    CAJFCJ010000053.1    19349    19649  +           cds-CAD5126493.1
       14  |    CAJFCJ010000053.1    38860    39160  +           cds-CAD5126495.1
  PyRanges with 5 rows, 5 columns, and 1 index columns.
  Contains 2 chromosomes and 2 strands.

To detect cases of out-of-bounds on the right side, :func:`genome_bounds <pyranges.PyRanges.genome_bounds>`
needs to know chromosome sizes.
Various input types are accepted for the ``chromsizes`` argument; we used a ``pyfaidx.Fasta``
object, which derives it from a fasta file.

You see below that some intervals were gone out-of-bounds on the right side, and have been corrected:

  >>> diff_end = cor_prom.End != prom.End
  >>> prom[diff_end]
    index  |    Chromosome           Start      End  Strand      ID
    int64  |    category             int64    int64  category    object
  -------  ---  -----------------  -------  -------  ----------  ----------------
        4  |    CAJFCJ010000025.1     3153     3453  -           cds-CAD5125114.1
  PyRanges with 1 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.

  >>> cor_prom[diff_end]
    index  |    Chromosome           Start      End  Strand      ID
    int64  |    category             int64    int64  category    object
  -------  ---  -----------------  -------  -------  ----------  ----------------
        4  |    CAJFCJ010000025.1     3153     3418  -           cds-CAD5125114.1
  PyRanges with 1 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.


Detecting overlaps among intervals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pyranges offers many efficient methods to detect overlaps, such as
:func:`overlap <pyranges.PyRanges.overlap>`.
This method returns the rows in self that overlap with another PyRanges object.

Let's see if any of the promoter regions overlap other CDSs:

  >>> cor_prom.overlap(cds)
    index  |    Chromosome           Start      End  Strand      ID
    int64  |    category             int64    int64  category    object
  -------  ---  -----------------  -------  -------  ----------  ----------------
        0  |    CAJFCJ010000025.1     2755     3055  -           cds-CAD5125115.1
  PyRanges with 1 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.

As many PyRanges methods, the Strand (if present) is taken into account in the comparison, so that
the overlap bewteen intervals is reported only if they are on the same strand.
Argument ``strand_behavior`` is available in many functions to control how strand is handled in overlap comparisons.
(see :func:`overlap <pyranges.PyRanges.overlap>`).

Above, we obtained the promoter region that overlaps another CDS, but we don't know what CDS it is.
Function :func:`join_ranges <pyranges.PyRanges.join_ranges>` will find overlaps and combine the columns
of the overlapping intervals, similar to a SQL join operation:

  >>> j = cor_prom.join_ranges(cds)
  >>> j
    index  |    Chromosome           Start      End  Strand      ID                  Start_b    End_b  ID_b
    int64  |    category             int64    int64  category    object                int64    int64  object
  -------  ---  -----------------  -------  -------  ----------  ----------------  ---------  -------  ----------------
        0  |    CAJFCJ010000025.1     2755     3055  -           cds-CAD5125115.1       2753     2851  cds-CAD5125114.1
  PyRanges with 1 rows, 8 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.

The object ``j`` contains the columns of both objects, with the suffix "_b" to distinguish the second one (``cds``).
It may be a bit too wide for our taste. Let's just look at a few columns to understand the overlap:

  >>> j[['ID', 'Start', 'End', 'ID_b', 'Start_b', 'End_b']]
                   ID  Start   End              ID_b  Start_b  End_b
  0  cds-CAD5125115.1   2755  3055  cds-CAD5125114.1     2753   2851

Above, we used a pandas syntax to select columns. Because the returned object does not have all genomic location
columns, it is a pandas DataFrame.

Let's get the intersection between the overlapping intervals, using function
:func:`intersect <pyranges.PyRanges.intersect>`:

  >>> prom_in_cds = cor_prom.intersect(cds)
  >>> prom_in_cds
    index  |    Chromosome           Start      End  Strand      ID
    int64  |    category             int64    int64  category    object
  -------  ---  -----------------  -------  -------  ----------  ----------------
        0  |    CAJFCJ010000025.1     2755     2851  -           cds-CAD5125115.1
  PyRanges with 1 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.

Let's go back to the ``cds`` object and see if any of its intervals overlap each other.
We can use :func:`cluster <pyranges.PyRanges.cluster>`. This will assign each interval to a cluster,
identified by an integer. The intervals that overlap each other will be assigned to the same cluster.

  >>> clu_cds = cds.cluster()
  >>> clu_cds
  index    |    Chromosome         Start    End      Strand      ID                Cluster
  int64    |    category           int64    int64    category    object            int64
  -------  ---  -----------------  -------  -------  ----------  ----------------  ---------
  4        |    CAJFCJ010000053.1  4882     5263     -           cds-CAD5126491.1  0
  11       |    CAJFCJ010000053.1  10732    10958    +           cds-CAD5126492.1  1
  12       |    CAJFCJ010000053.1  11028    11169    +           cds-CAD5126492.1  2
  13       |    CAJFCJ010000053.1  11227    11400    +           cds-CAD5126492.1  3
  ...      |    ...                ...      ...      ...         ...               ...
  146      |    CAJFCJ010000025.1  2753     2851     -           cds-CAD5125114.1  42
  147      |    CAJFCJ010000025.1  2593     2693     -           cds-CAD5125114.1  43
  148      |    CAJFCJ010000025.1  2354     2537     -           cds-CAD5125114.1  44
  149      |    CAJFCJ010000025.1  2174     2294     -           cds-CAD5125114.1  45
  PyRanges with 56 rows, 6 columns, and 1 index columns.
  Contains 3 chromosomes and 2 strands.

Let's get the clusters that have more than one interval in them, using pandas ``value_counts``:

  >>> c = clu_cds.Cluster.value_counts()
  >>> multi_clusters = c[ c > 1 ].index
  >>> multi_clu_cds = clu_cds[ clu_cds.Cluster.isin(multi_clusters) ]
  >>> multi_clu_cds
  index    |    Chromosome         Start    End      Strand      ID                Cluster
  int64    |    category           int64    int64    category    object            int64
  -------  ---  -----------------  -------  -------  ----------  ----------------  ---------
  110      |    CAJFCJ010000097.1  51865    52382    +           cds-CAD5126878.1  38
  111      |    CAJFCJ010000097.1  52446    52826    +           cds-CAD5126878.1  39
  112      |    CAJFCJ010000097.1  52903    53027    +           cds-CAD5126878.1  40
  113      |    CAJFCJ010000097.1  53339    53404    +           cds-CAD5126878.1  41
  ...      |    ...                ...      ...      ...         ...               ...
  146      |    CAJFCJ010000025.1  2753     2851     -           cds-CAD5125114.1  42
  147      |    CAJFCJ010000025.1  2593     2693     -           cds-CAD5125114.1  43
  148      |    CAJFCJ010000025.1  2354     2537     -           cds-CAD5125114.1  44
  149      |    CAJFCJ010000025.1  2174     2294     -           cds-CAD5125114.1  45
  PyRanges with 17 rows, 6 columns, and 1 index columns.
  Contains 2 chromosomes and 2 strands.

Sorting intervals
~~~~~~~~~~~~~~~~~
Above, it is not apparent that there are overlaps among the intervals in the object ``multi_clu_cds``. This is due to
the order of rows. We could sort row using pandas ``sort_values``, but PyRanges offers something
better: the method :func:`sort_ranges <pyranges.PyRanges.sort_ranges>` sorts by chromosome, strand, then by
coordinates. By default, intervals are sorted 5' to 3', meaning that intervals on the positive strand are sorted
from left-most to right-most, while intervals on the negative strand are sorted in the opposite direction.

  >>> multi_clu_cds.sort_ranges()
  index    |    Chromosome         Start    End      Strand      ID                Cluster
  int64    |    category           int64    int64    category    object            int64
  -------  ---  -----------------  -------  -------  ----------  ----------------  ---------
  146      |    CAJFCJ010000025.1  2753     2851     -           cds-CAD5125114.1  42
  135      |    CAJFCJ010000025.1  2753     2755     -           cds-CAD5125115.1  42
  136      |    CAJFCJ010000025.1  2593     2693     -           cds-CAD5125115.1  43
  147      |    CAJFCJ010000025.1  2593     2693     -           cds-CAD5125114.1  43
  ...      |    ...                ...      ...      ...         ...               ...
  112      |    CAJFCJ010000097.1  52903    53027    +           cds-CAD5126878.1  40
  123      |    CAJFCJ010000097.1  52903    53027    +           cds-CAD5126877.1  40
  113      |    CAJFCJ010000097.1  53339    53404    +           cds-CAD5126878.1  41
  124      |    CAJFCJ010000097.1  53339    53404    +           cds-CAD5126877.1  41
  PyRanges with 17 rows, 6 columns, and 1 index columns.
  Contains 2 chromosomes and 2 strands.


:func:`sort_ranges <pyranges.PyRanges.sort_ranges>` offers many options to customize the sorting.
For example, let's add a columns with the lengths of each interval.
Thus, sort by chromosome, strand, length (in descending order), then interval coordinates:

  >>> multi_clu_cds['Length'] = multi_clu_cds.lengths()
  >>> multi_clu_cds.sort_ranges(by='Length', sort_descending='Length')
  index    |    Chromosome         Start    End      Strand      ID                Cluster    Length
  int64    |    category           int64    int64    category    object            int64      int64
  -------  ---  -----------------  -------  -------  ----------  ----------------  ---------  --------
  137      |    CAJFCJ010000025.1  2354     2537     -           cds-CAD5125115.1  44         183
  148      |    CAJFCJ010000025.1  2354     2537     -           cds-CAD5125114.1  44         183
  138      |    CAJFCJ010000025.1  2174     2294     -           cds-CAD5125115.1  45         120
  149      |    CAJFCJ010000025.1  2174     2294     -           cds-CAD5125114.1  45         120
  ...      |    ...                ...      ...      ...         ...               ...        ...
  123      |    CAJFCJ010000097.1  52903    53027    +           cds-CAD5126877.1  40         124
  121      |    CAJFCJ010000097.1  52261    52382    +           cds-CAD5126877.1  38         121
  113      |    CAJFCJ010000097.1  53339    53404    +           cds-CAD5126878.1  41         65
  124      |    CAJFCJ010000097.1  53339    53404    +           cds-CAD5126877.1  41         65
  PyRanges with 17 rows, 7 columns, and 1 index columns.
  Contains 2 chromosomes and 2 strands.

Other overlap-based operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Say that we are interested in intergenic regions in chromosome ``CAJFCJ010000097.1``.
In any genome annotation, the annotation rows "exon" define transcript coordinates. Let's fetch them from
the annotation ``ann``:

  >>> exons = ann[ ann.Feature == 'exon' ].loci['CAJFCJ010000097.1']
  >>> exons
  index    |    Chromosome         Start    End      Strand      Feature     Parent                 ...
  int64    |    category           int64    int64    category    category    object                 ...
  -------  ---  -----------------  -------  -------  ----------  ----------  ---------------------  -----
  86       |    CAJFCJ010000097.1  2248     3308     +           exon        rna-DGYR_LOCUS14091    ...
  90       |    CAJFCJ010000097.1  6505     6600     -           exon        rna-DGYR_LOCUS14092    ...
  91       |    CAJFCJ010000097.1  6082     6450     -           exon        rna-DGYR_LOCUS14092    ...
  92       |    CAJFCJ010000097.1  5579     6029     -           exon        rna-DGYR_LOCUS14092    ...
  ...      |    ...                ...      ...      ...         ...         ...                    ...
  116      |    CAJFCJ010000097.1  52261    52382    +           exon        rna-DGYR_LOCUS14095-2  ...
  117      |    CAJFCJ010000097.1  52446    52826    +           exon        rna-DGYR_LOCUS14095-2  ...
  118      |    CAJFCJ010000097.1  52903    53027    +           exon        rna-DGYR_LOCUS14095-2  ...
  119      |    CAJFCJ010000097.1  53339    53404    +           exon        rna-DGYR_LOCUS14095-2  ...
  PyRanges with 15 rows, 8 columns, and 1 index columns. (2 columns not shown: "ID", "midpoint").
  Contains 1 chromosomes and 2 strands.

Let's define the boundaries of each mRNA, e.g. the left and right limits of its exons. While this may be readily
available in the genome annotation, let's use PyRanges to calculate them, using
:func:`boundaries <pyranges.PyRanges.boundaries>`:

  >>> mRNA_bounds = exons.boundaries(transcript_id='Parent')
  >>> mRNA_bounds
    index  |    Chromosome           Start      End  Strand      Parent
    int64  |    category             int64    int64  category    object
  -------  ---  -----------------  -------  -------  ----------  ---------------------
        0  |    CAJFCJ010000097.1     2248     3308  +           rna-DGYR_LOCUS14091
        1  |    CAJFCJ010000097.1    16697    17634  +           rna-DGYR_LOCUS14093
        2  |    CAJFCJ010000097.1    51864    53404  +           rna-DGYR_LOCUS14095
        3  |    CAJFCJ010000097.1    51864    53404  +           rna-DGYR_LOCUS14095-2
        4  |    CAJFCJ010000097.1     5579     6600  -           rna-DGYR_LOCUS14092
        5  |    CAJFCJ010000097.1    31876    32195  -           rna-DGYR_LOCUS14094
  PyRanges with 6 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 2 strands.

To get the intergenic regions, let's define the maximum and minimum coordinates of any mRNA in this region,
using :func:`boundaries <pyranges.PyRanges.boundaries>` again without ``transcript_id``. Because we want our result to
not depend on strand, we remove it using :func:`remove_strand <pyranges.PyRanges.remove_strand>`:

  >>> all_mRNA_bounds = mRNA_bounds.remove_strand().boundaries()
  >>> all_mRNA_bounds
    index  |    Chromosome           Start      End
    int64  |    category             int64    int64
  -------  ---  -----------------  -------  -------
        0  |    CAJFCJ010000097.1     2248    53404
  PyRanges with 1 rows, 3 columns, and 1 index columns.
  Contains 1 chromosomes.

Now we can get the intergenic regions using :func:`subtract_ranges <pyranges.PyRanges.subtract_ranges>`:

  >>> intergenic = all_mRNA_bounds.subtract_ranges(mRNA_bounds)
  >>> intergenic
    index  |    Chromosome           Start      End
    int64  |    category             int64    int64
  -------  ---  -----------------  -------  -------
        0  |    CAJFCJ010000097.1     3308     5579
        0  |    CAJFCJ010000097.1     6600    16697
        0  |    CAJFCJ010000097.1    17634    31876
        0  |    CAJFCJ010000097.1    32195    51864
  PyRanges with 4 rows, 3 columns, and 1 index columns (with 3 index duplicates).
  Contains 1 chromosomes.

Note that pyranges indicates that the object has duplicate indices, because all come from the same row in
``all_mRNA_bounds``, broken into subintervals by the subtraction operation.
We can use pandas ``reset_index`` to remedy:

  >>> intergenic = intergenic.reset_index(drop=True)
  >>> intergenic
    index  |    Chromosome           Start      End
    int64  |    category             int64    int64
  -------  ---  -----------------  -------  -------
        0  |    CAJFCJ010000097.1     3308     5579
        1  |    CAJFCJ010000097.1     6600    16697
        2  |    CAJFCJ010000097.1    17634    31876
        3  |    CAJFCJ010000097.1    32195    51864
  PyRanges with 4 rows, 3 columns, and 1 index columns.
  Contains 1 chromosomes.


Counting overlaps
~~~~~~~~~~~~~~~~~

Often, one wants to count the number of overlaps between two PyRanges objects, e.g. to count reads in specific regions.
Here, let's count the number of CDS intervals that overlap our previously computed objects ``intergenic``
and  ``all_mRNA_bounds``, using  method :func:`count_overlaps <pyranges.PyRanges.count_overlaps>` :

  >>> intergenic.count_overlaps(cds)
    index  |    Chromosome           Start      End    NumberOverlaps
    int64  |    category             int64    int64             int64
  -------  ---  -----------------  -------  -------  ----------------
        0  |    CAJFCJ010000097.1     3308     5579                 0
        1  |    CAJFCJ010000097.1     6600    16697                 0
        2  |    CAJFCJ010000097.1    17634    31876                 0
        3  |    CAJFCJ010000097.1    32195    51864                 0
  PyRanges with 4 rows, 4 columns, and 1 index columns.
  Contains 1 chromosomes.

  >>> all_mRNA_bounds.count_overlaps(cds)
    index  |    Chromosome           Start      End    NumberOverlaps
    int64  |    category             int64    int64             int64
  -------  ---  -----------------  -------  -------  ----------------
        0  |    CAJFCJ010000097.1     2248    53404                15
  PyRanges with 1 rows, 4 columns, and 1 index columns.
  Contains 1 chromosomes.

As expected, there's no CDS overlapping the intergenic regions, while the other object reports 15. Yet,
the CDS intervals may be redundant: different splicing isoforms may have some identical exons:

  >>> example = cds.loci['CAJFCJ010000097.1', '+', 51000:54000].sort_ranges()
  >>> example
  index    |    Chromosome         Start    End      Strand      ID
  int64    |    category           int64    int64    category    object
  -------  ---  -----------------  -------  -------  ----------  ----------------
  120      |    CAJFCJ010000097.1  51865    52201    +           cds-CAD5126877.1
  110      |    CAJFCJ010000097.1  51865    52382    +           cds-CAD5126878.1
  121      |    CAJFCJ010000097.1  52261    52382    +           cds-CAD5126877.1
  111      |    CAJFCJ010000097.1  52446    52826    +           cds-CAD5126878.1
  ...      |    ...                ...      ...      ...         ...
  112      |    CAJFCJ010000097.1  52903    53027    +           cds-CAD5126878.1
  123      |    CAJFCJ010000097.1  52903    53027    +           cds-CAD5126877.1
  113      |    CAJFCJ010000097.1  53339    53404    +           cds-CAD5126878.1
  124      |    CAJFCJ010000097.1  53339    53404    +           cds-CAD5126877.1
  PyRanges with 9 rows, 5 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.


If we want to calculate the intervals that are annotated as CDS in any of the isoforms, we can use
method :func:`merge_overlaps <pyranges.PyRanges.merge_overlaps>` :

  >>> example.merge_overlaps()
    index  |    Chromosome           Start      End  Strand
    int64  |    object               int64    int64  object
  -------  ---  -----------------  -------  -------  --------
        0  |    CAJFCJ010000097.1    51865    52382  +
        1  |    CAJFCJ010000097.1    52446    52826  +
        2  |    CAJFCJ010000097.1    52903    53027  +
        3  |    CAJFCJ010000097.1    53339    53404  +
  PyRanges with 4 rows, 4 columns, and 1 index columns.
  Contains 1 chromosomes and 1 strands.

Various methods are available to obtain non-overlapping intervals, depending on the desired output. See
:func:`split <pyranges.PyRanges.split>`, :func:`max_disjoint <pyranges.PyRanges.max_disjoint>`.

Finally, let's count how many non-redundant CDS intervals overlap our target region:

  >>> all_mRNA_bounds.count_overlaps(cds.merge_overlaps())
    index  |    Chromosome           Start      End    NumberOverlaps
    int64  |    category             int64    int64             int64
  -------  ---  -----------------  -------  -------  ----------------
        0  |    CAJFCJ010000097.1     2248    53404                10
  PyRanges with 1 rows, 4 columns, and 1 index columns.
  Contains 1 chromosomes.


This concludes our tutorial. The next pages will delve into pyranges functionalities grouped by topic.
