Mapping between coordinate systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. contents::
   :local:
   :depth: 2

Mapping coordinate systems: cheatsheet
--------------------------------------

.. image:: https://raw.githubusercontent.com/pyranges/pyranges_plot/for_pyranges1_1/examples/cheatsheet_mapping.png
   :alt: PyRanges cheatsheet
   :target: https://raw.githubusercontent.com/pyranges/pyranges_plot/for_pyranges1_1/examples/cheatsheet_mapping.png


:func:`overlap <pyranges.PyRanges.overlap>`.

Nested coordinate systems
-------------------------

PyRanges may represent nested coordinate systems, where the intervals in one coordinate system
are relative to the intervals in another coordinate system. This is useful for representing hierarchical data.

Let's make an intuitive case. On the one hand, we have some global ranges, that is
genomic intervals wherein coordinates refer to positions along the genome:

  >>> import pyranges as pr
  >>> gr = pr.example_data.ncbi_gff
  >>> gre = (gr[gr.Feature=='exon']).get_with_loc_columns('ID')
  >>> gre['ID']= gre['ID'].str.replace('exon-', '', regex=False) # remove the 'exon-' prefix for clarity
  >>> gre
  index    |    Chromosome         Start    End      Strand      ID
  int64    |    category           int64    int64    category    object
  -------  ---  -----------------  -------  -------  ----------  -------------------
  3        |    CAJFCJ010000053.1  4882     5264     -           DGYR_LOCUS13733-1
  7        |    CAJFCJ010000053.1  10474    10958    +           DGYR_LOCUS13734-1
  8        |    CAJFCJ010000053.1  11028    11169    +           DGYR_LOCUS13734-2
  9        |    CAJFCJ010000053.1  11227    11400    +           DGYR_LOCUS13734-3
  ...      |    ...                ...      ...      ...         ...
  141      |    CAJFCJ010000025.1  2753     2851     -           DGYR_LOCUS12552-2-2
  142      |    CAJFCJ010000025.1  2593     2693     -           DGYR_LOCUS12552-2-3
  143      |    CAJFCJ010000025.1  2354     2537     -           DGYR_LOCUS12552-2-4
  144      |    CAJFCJ010000025.1  1909     2294     -           DGYR_LOCUS12552-2-5
  PyRanges with 57 rows, 5 columns, and 1 index columns.
  Contains 3 chromosomes and 2 strands.

On the other hand, we have local data, i.e. mapped to the transcript sequences that are annotated in the global ranges.
Below, we load some example data that we obtained by running the Rfam database of RNA motifs against the
transcript sequences in the ``gre`` object, using the software Infernal:

  >>> rh = pr.example_data.rfam_hits
  >>> rh = rh[ ['target_name', 'seq_from', 'seq_to', 'strand', 'query_name', 'mdl_from', 'mdl_to'] ]
  >>> rh.head()
             target_name  seq_from  seq_to strand query_name  mdl_from  mdl_to
  0  DGYR_LOCUS12552-2-3        51      95      +       GAIT         1      71
  1    DGYR_LOCUS12552-3        51      95      +       GAIT         1      71
  2    DGYR_LOCUS13738-3        87     124      +    REN-SRE         1      37
  3    DGYR_LOCUS14091-1       137      43      -  IFN_gamma         1     169
  4    DGYR_LOCUS14091-1       547     616      +     snoZ30         1      97

Above, the ``target_name`` column contains the transcript IDs, which are
the same as the ``ID`` column in the ``gre`` object. Then we have coordinates that define the alignment between
portions of those transcripts and one of the RNA motifs (query models) in the Rfam database.

Let's convert the ``rh`` object to a PyRanges object, taking care of switching from 1-based to 0-based coordinates:

  >>> rh = rh.rename(columns={'target_name':'Chromosome', 'seq_from':'Start', 'seq_to':'End', 'strand':'Strand'})
  >>> rh = pr.PyRanges(rh)
  >>> rh['Start'] -= 1  # convert to 0-based coordinates
  >>> rh
  index    |    Chromosome           Start    End      Strand    query_name    mdl_from    mdl_to
  int64    |    object               int64    int64    object    object        int64       int64
  -------  ---  -------------------  -------  -------  --------  ------------  ----------  --------
  0        |    DGYR_LOCUS12552-2-3  50       95       +         GAIT          1           71
  1        |    DGYR_LOCUS12552-3    50       95       +         GAIT          1           71
  2        |    DGYR_LOCUS13738-3    86       124      +         REN-SRE       1           37
  3        |    DGYR_LOCUS14091-1    136      43       -         IFN_gamma     1           169
  ...      |    ...                  ...      ...      ...       ...           ...         ...
  30       |    DGYR_LOCUS13737-1    549      600      +         mir-3047      1           61
  31       |    DGYR_LOCUS13737-1    543      605      +         mir-3156      1           77
  32       |    DGYR_LOCUS13737-1    529      612      +         MIR1523       1           92
  33       |    DGYR_LOCUS13737-1    549      600      +         MIR8001       1           67
  PyRanges with 34 rows, 7 columns, and 1 index columns.
  Contains 18 chromosomes and 2 strands.
  Invalid ranges:
    * 14 intervals are empty or negative length (end <= start). See indexes: 3, 5, 6, ...

Infernal, like other programs, reports negative stranded hits with Start and End reversed. Let's fix that:

  >>> rh.loc[ rh.Strand == '-', ['Start', 'End'] ] = rh.loc[ rh.Strand == '-', ['End', 'Start'] ].values
  >>> rh
  index    |    Chromosome           Start    End      Strand    query_name    mdl_from    mdl_to
  int64    |    object               int64    int64    object    object        int64       int64
  -------  ---  -------------------  -------  -------  --------  ------------  ----------  --------
  0        |    DGYR_LOCUS12552-2-3  50       95       +         GAIT          1           71
  1        |    DGYR_LOCUS12552-3    50       95       +         GAIT          1           71
  2        |    DGYR_LOCUS13738-3    86       124      +         REN-SRE       1           37
  3        |    DGYR_LOCUS14091-1    43       136      -         IFN_gamma     1           169
  ...      |    ...                  ...      ...      ...       ...           ...         ...
  30       |    DGYR_LOCUS13737-1    549      600      +         mir-3047      1           61
  31       |    DGYR_LOCUS13737-1    543      605      +         mir-3156      1           77
  32       |    DGYR_LOCUS13737-1    529      612      +         MIR1523       1           92
  33       |    DGYR_LOCUS13737-1    549      600      +         MIR8001       1           67
  PyRanges with 34 rows, 7 columns, and 1 index columns.
  Contains 18 chromosomes and 2 strands.

Now we have the ``gre`` and ``rh`` objects, that represent global and local coordinate systems, respectively.
Let's check that all ``Chromosome`` values in ``rh`` matches the ``ID`` column in ``gre``:
  >>> bool( rh.Chromosome.isin(gre.ID).all() )
  True

Mapping from local to global ranges
-----------------------------------

Next, we want to take the Rfam hits in ``rh``, which are relative (local) to the transcript sequences,
and remap them to the genome (global) coordinates. To do so, we make use of the information in the ``gre`` object,
which defines the coordinates of each transcript, often split in exons, relative to the genome.

For this operation, we use the :func:`map_to_global <pyranges.PyRanges.map_to_global>` method.
Besides the two PyRanges objects, we also need to specify the columns in the global range
that contains the identifier used as the Chromosome in the local range:

  >>> rhg = rh.map_to_global(gre, global_on='ID')
  >>> rhg
  index    |    Chromosome         Start    End      Strand    query_name    mdl_from    mdl_to
  int64    |    category           int64    int64    object    object        int64       int64
  -------  ---  -----------------  -------  -------  --------  ------------  ----------  --------
  0        |    CAJFCJ010000025.1  2598     2643     -         GAIT          1           71
  1        |    CAJFCJ010000025.1  2598     2643     -         GAIT          1           71
  2        |    CAJFCJ010000053.1  77544    77582    +         REN-SRE       1           37
  3        |    CAJFCJ010000097.1  2291     2384     -         IFN_gamma     1           169
  ...      |    ...                ...      ...      ...       ...           ...         ...
  11       |    CAJFCJ010000097.1  5875     5921     -         snR77         17          61
  12       |    CAJFCJ010000097.1  5700     5777     -         snR50         1           89
  21       |    CAJFCJ010000097.1  51976    52037    -         MESTIT1_1     60          124
  22       |    CAJFCJ010000097.1  51976    52037    -         MESTIT1_1     60          124
  PyRanges with 9 rows, 7 columns, and 1 index columns.
  Contains 3 chromosomes and 2 strands.

// then mention: keep columns
// modify by adding an example that maps to two exons, add to example data

// placeholders for later:



  >>> grc = (gr[gr.Feature=='CDS']).get_with_loc_columns(['ID', 'protein_id'])
  >>> grc
  index    |    Chromosome         Start    End      Strand      ID                protein_id
  int64    |    category           int64    int64    category    object            object
  -------  ---  -----------------  -------  -------  ----------  ----------------  ------------
  4        |    CAJFCJ010000053.1  4882     5263     -           cds-CAD5126491.1  CAD5126491.1
  11       |    CAJFCJ010000053.1  10732    10958    +           cds-CAD5126492.1  CAD5126492.1
  12       |    CAJFCJ010000053.1  11028    11169    +           cds-CAD5126492.1  CAD5126492.1
  13       |    CAJFCJ010000053.1  11227    11400    +           cds-CAD5126492.1  CAD5126492.1
  ...      |    ...                ...      ...      ...         ...               ...
  146      |    CAJFCJ010000025.1  2753     2851     -           cds-CAD5125114.1  CAD5125114.1
  147      |    CAJFCJ010000025.1  2593     2693     -           cds-CAD5125114.1  CAD5125114.1
  148      |    CAJFCJ010000025.1  2354     2537     -           cds-CAD5125114.1  CAD5125114.1
  149      |    CAJFCJ010000025.1  2174     2294     -           cds-CAD5125114.1  CAD5125114.1
  PyRanges with 56 rows, 6 columns, and 1 index columns.
  Contains 3 chromosomes and 2 strands.


  >>> genome_file = pr.example_data.files['ncbi.fasta']
  >>> mrna_seq = gre.get_sequence(genome_file, group_by='ID').str.upper()
  >>> cds_seq = gr.get_sequence(genome_file, group_by='ID').str.upper()
  >>> pep_seq = pr.seqs.translate(cds_seq)
