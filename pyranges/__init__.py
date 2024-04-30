# must import here every module/object X that will be available as pyranges.X
# For each 'extension' e.g. named foobar, you must have:
# a ext/foobar.py file in pyranges/ which contains the actual code for methods etc
# a foobar.py file in pyranges/ which just imports the minimal objects to be exposed

from pyranges import orfs, seqs, stats  # noqa: F401
from pyranges.core.example_data import example_data  # noqa: F401
from pyranges.core.multioverlap import count_overlaps  # noqa: F401
from pyranges.core.options import option_manager as options  # noqa: F401
from pyranges.core.pyranges_main import PyRanges  # noqa: F401
from pyranges.core.random import random  # noqa: F401
from pyranges.core.version import __version__  # noqa: F401
from pyranges.methods.concat import concat  # noqa: F401
from pyranges.methods.tile_genome import tile_genome  # noqa: F401
from pyranges.range_frame.range_frame import RangeFrame  # noqa: F401
from pyranges.readers import from_string, read_bam, read_bed, read_bigwig, read_gff3, read_gtf  # noqa: F401

read_gff = read_gtf
