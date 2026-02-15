.. pyranges1 documentation master file, created by
   sphinx-quickstart.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


The pyranges1 documentation
==========================
Pyranges is a Python library specifically designed for efficient and intuitive manipulation of genomics data,
particularly genomic intervals (like genes, genomic features, or reads).
The library is optimized for fast querying and manipulation of genomic annotations.

Pyranges is open source, and hosted at github: https://github.com/pyranges/pyranges_1.x/

Pyranges is developed by Endre Bakken Stovner and by
`Marco Mariotti's lab <https://www.mariottigenomicslab.com/>`_.

**This documentation refers to the version 1 of pyranges**, which introduces a new API and a new data structure
compared to version 0, still the 'default', soon to be deprecated. If you are a pyranges1 v0 user,
check :doc:`this guide to migrate to v1 <./migration_guide>`.


While we recommend using pyranges1 as Python library, alternatively, the
:doc:`pyranger command-line tool <./pyranger_cli>` allows you to access the same functionalities from the command line,
without writing any Python code.

Citation
~~~~~~~~

Stovner EB, Ticó M, Muñoz del Campo E, Pallarès-Albanell J, Chawla K, Sætrom P, Mariotti M (2025) Pyranges v1: a Python framework for ultrafast sequence interval operations.
*bioRxiv* 2025.12.11.693639; doi: https://doi.org/10.64898/2025.12.11.693639

Stovner EB, Sætrom P (2020) PyRanges: efficient comparison of genomic intervals in Python. *Bioinformatics 36(3):918-919*  http://dx.doi.org/10.1093/bioinformatics/btz615


Documentation outline
~~~~~~~~~~~~~~~~~~~~~

#. 🛠️  :doc:`Installation instructions <./installation>`
#. 📖  :doc:`The tutorial <./tutorial>`,  recommended for all new users
#. 🧰  :doc:`The how-to pages <./how_to_pages>`, where functionalities are grouped by topic
#. 📑  :doc:`The API reference <./api_reference>`, where all methods are explained in detail
#. 🚀  :doc:`The cheatsheet <./cheatsheet>`, a quick reference for the most common operations
#. 💻  :doc:`The pyranger command-line tool <./pyranger_cli>`, for those who prefer the command line
#. 🧑‍💻 :doc:`The developer guide <./developer_guide>`, to follow in order to contribute to PyRanges
#. 🔄  :doc:`The guide to migrate to v1 <./migration_guide>`, for existing users of PyRanges v0


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Contents:

   installation
   tutorial
   how_to_pages
   api_reference
   cheatsheet
   pyranger_cli
   developer_guide
   migration_guide


Indices
=======

* :ref:`genindex`
* :ref:`modindex`
