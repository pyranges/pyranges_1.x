import logging
from collections import OrderedDict
from typing import TYPE_CHECKING

import pandas as pd
from tabulate import tabulate

if TYPE_CHECKING:
    from pyranges import PyRanges


logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)
LOGGER.Formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s %(message)s")
LOGGER.setLevel(logging.INFO)

def _summary(
    self: "PyRanges",
    *,
    return_df: bool = False,
) -> pd.DataFrame | None:
    lengths = {}
    total_lengths = {}
    lengths["pyrange"] = self.lengths()
    total_lengths["pyrange"] = [self.length]

    if self.strand_values_valid:
        c = self.merge_overlaps(strand=True)
        lengths["coverage_forward"] = c.loci["+"].lengths()
        lengths["coverage_reverse"] = c.loci["-"].lengths()
        total_lengths["coverage_forward"] = [c.loci["+"].length]
        total_lengths["coverage_reverse"] = [c.loci["-"].length]
    else:
        c = self

    c = c.merge_overlaps(strand=False)
    lengths["coverage_unstranded"] = c.lengths()
    total_lengths["coverage_unstranded"] = [c.length]

    summaries = OrderedDict()

    # statistics for lengths
    for summary, s in lengths.items():
        summaries[summary] = s.describe()

    summary = pd.concat(summaries.values(), axis=1)
    summary.columns = list(summaries)

    df = pd.DataFrame.from_dict(total_lengths)
    df.index = ["sum"]
    summary = pd.concat([summary, df])

    if return_df:
        return summary
    str_repr = tabulate(summary, headers=summary.columns)
    print(str_repr)  # noqa: T201
    return None
