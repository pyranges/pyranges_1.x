# pyranges

## Introduction

Pyranges is a Python library with a Rust backend for efficient and intuitive manipulation of genomics data,
particularly genomic intervals (like genes, genomic features, or reads).
The library is optimized for fast querying and manipulation of genomic annotations.
It enables intuitive and highly efficient pipelines for genomic analysis.

*"Finally ... This was what Python badly needed for years."* - Heng Li

## Version 1.x
This is version 1.x of pyranges. It is a complete rewrite of the original pyranges library, 
that will replace the "default" (version 0) at the end in 2025. If you are a v0 user, check the migration guide 
in the documentation.

## Documentation

The pyranges documentation, including installation instructions, API, tutorial, and how-to-pages, is 
available at https://pyranges1.readthedocs.io/

## Recent Changelog

```
# 1.1.9 (26.01.26)
- pandas dependency bound to v2. This is in response to pandas 3.0.0 being released, breaking our doctests.

# 1.1.8 (30.12.25)
- to_gtf and to_gff3: fix bug where 'phase' (gtf) and 'frame' (gff3) are erroneously added to attributes field

# v1.1.7 (16.12.25)
- window_ranges: fix sort order issue in  when using by (#98 and #105)
- window_ranges: added argument add_window_id, updated documentation
```

## Install

Pyranges 1.x requires python ≥3.12. Minimal installation: 

```bash
pip install pyranges1
```

Installation including all optional dependencies:

```bash
pip install pyranges1[all]
```

Details at https://pyranges1.readthedocs.io/en/latest/installation.html


## Features

  - fast
  - memory-efficient
  - featureful
  - pythonic/pandastic

## Paper/Cite

For v1:

Stovner EB, Ticó M, Muñoz del Campo E, Pallarès-Albanell J, Chawla K, Sætrom P, Mariotti M (2025) Pyranges v1: a Python framework for ultrafast sequence interval operations.
*bioRxiv* 2025.12.11.693639; doi: https://doi.org/10.64898/2025.12.11.693639


For v0:

Stovner EB, Sætrom P (2020) PyRanges: efficient comparison of genomic intervals in Python. 
*Bioinformatics 36(3):918-919*  http://dx.doi.org/10.1093/bioinformatics/btz615

## Supporting pyranges

  - most importantly, cite pyranges if you use it. It is the main metric funding sources care about.
  - use pyranges in Stack Overflow/biostars questions and answers
  - star the repo (possibly important for github visibility and as a proxy for project popularity)

## Asking for help

If you encounter bugs, or the documentation is not enough a cannot accomplish a specific task of interest, 
or if you'd like new features implemented, open an Issue at github: https://github.com/pyranges/pyranges/issues

## Contributing to pyranges

Pyranges accepts code contributions in form of pull request. 
For details, visit https://pyranges1.readthedocs.io/developer_guide.html

## Cheatsheet
![cheatsheet](https://raw.githubusercontent.com/pyranges/pyrangeyes/for_pyranges1_1/images/pyranges_cheatsheet.png)
(The cheatsheet above was created with pyrangeyes, a companion graphical library:  https://pyrangeyes.readthedocs.io/)




