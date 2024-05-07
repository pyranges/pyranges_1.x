Migration guide from v0 to v1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pyranges v1.x introduced major changes in data structure and interface.
This guide is intended to help you migrate your code from v0.x to v1.x.

Useful links:

 * `Documentation for v0 <https://pyranges.readthedocs.io/en/v0/>`_
 * `Github repository for v0 <https://github.com/pyranges/pyranges/tree/v0>`_

In v0, PyRanges objects were implemented as a dictionary of DataFrames, one for each chromosome or for each
(chromosome, strand) pair if strand was present and valid. This implied that intervals on different chromosomes had
no inherent order, and operations on each chromosome were independent. In practical uses,
PyRanges had to be often converted to a single DataFrame to access the full functionality of pandas, using the
``df`` or ``as_df`` attributes, which was inefficient and cumbersome. Moreover, performance was not
very low on datasets comprising many "chromosomes" such as transcriptomic-based intervals.

In v1, PyRanges objects are implemented as a DataFrame subclass, which allows for more efficient operations and
direct access to pandas methods. This required a major interface change, since some methods and attributes existed
with the same name in both PyRanges and DataFrame, and the new implementation had to be consistent with the DataFrame.
Ultimately, we took advantage of the necessity of the change to redesign the interface to be more consistent and
maintainable in the future.

Below, we provide a cheatsheet to help you migrate your code from v0 to v1.
The most problematic aspects are the get/set item methods, since the v0 syntax is not compatible with dataframes.
See `here a discussion on the topic <https://github.com/pyranges/pyranges/discussions/357#discussioncomment-7274998>`_.

In the table below, ``pr`` refers to the pyranges module, and ``g`` to a PyRanges object. Most items are linked
to the corresponding documentation page.

.. csv-table:: Migration cheatsheet
   :file: migration_cheatsheet.tsv
   :delim: tab
   :header-rows: 1

