# must import here every module/object X that will be available as pyranges.X

# external (i.e. non-core) submodules must be imported here
# each will be available as pyranges.submodule_name

# note that there may be functions in external submodules that accept a PyRanges as first argument;
# these can be made available in pyranges.PyRanges.ext_submodule_name
# check all usages of 'StatsManager' to see how this is done

from pyranges import genomicfeatures  # noqa: F401
from pyranges.core.example_data_manager import example_data  # noqa: F401
from pyranges.core.multioverlap import count_overlaps  # noqa: F401
from pyranges.core.options import option_manager as options  # noqa: F401
from pyranges.core.pyranges_main import PyRanges  # noqa: F401
from pyranges.core.random import random  # noqa: F401
from pyranges.core.version import __version__  # noqa: F401
from pyranges.ext import stats  # noqa: F401
from pyranges.methods.concat import concat  # noqa: F401
from pyranges.range_frame.range_frame import RangeFrame  # noqa: F401
from pyranges.readers import from_string, read_bam, read_bed, read_bigwig, read_gff3, read_gtf  # NOQA: F401

read_gff = read_gtf
