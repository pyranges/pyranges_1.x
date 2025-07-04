
Installation
~~~~~~~~~~~~

Pyranges 1.x requires python ≥ 3.12. If necessary, create a new conda environment::

    conda create -yn pr python
    conda activate pr


The preferred way to install pyranges is via pip::

    pip install pyranges1

The command above will install a **minimal version** of pyranges.
Pyranges has several optional dependencies, required for certain functionalities.

To **install all optional dependencies**, use::

    pip install pyranges1[all]

Here you can see the optional dependencies grouped by functionality::

    # user add-ons: to fetch sequences, read BAM files, parallelize execution ...
    pip install pyranges1[add-ons]

    # command line: to use the pyranger command-line tool
    pip install pyranges1[cli]

    # development: for testing, linting, type checking, generating documentation
    pip install pyranges1[dev]

    # documentation: for building the documentation
    pip install pyranges1[docs]

To inspect the list of dependencies, check the pyproject.toml file in the repository.