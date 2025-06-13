The pyranges module
-------------------

The pyranges module exposes the class :class:`PyRanges <pyranges.PyRanges>` (omitted in this page)
as well as a number of functions for reading data from commonly used file formats.
You also have an :ref:`pyranges.options <pyranges_options>` interface to
configure how PyRanges objects are represented, and an :ref:`pyranges.example_data <pyranges_example_data>` object used in tests and documentation.


.. automodule:: pyranges
    :members:
    :imported-members:  # Ensure this is set to include imported members
    :exclude-members: PyRanges, RangeFrame

.. _pyranges_options:
pyranges.options
~~~~~~~~~~~~~~~~

The ``pyranges.options`` object is used to configure aspects of how PyRanges is represented.
Below are the methods available on this object.

.. automodule:: pyranges.options
   :members:
   :show-inheritance:
   :exclude-members: __init__

.. _pyranges_example_data:
pyranges.example_data
~~~~~~~~~~~~~~~~~~~~~
The ``pyranges.example_data`` object contains example data used in tests and documentation.
Printing it shows an overview of available data:

  >>> import pyranges as pr
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
  example_data.files              : A dict of basenames to file paths of available data.

Most of the data is in the form of PyRanges objects:

  >>> pr.example_data.chipseq
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
