from collections.abc import Iterable
from typing import TYPE_CHECKING, Any, Final, Literal, Protocol, TypeVar, get_args

if TYPE_CHECKING:
    import pandas as pd

    from pyranges import PyRanges, RangeFrame

RangeFrameType = TypeVar("RangeFrameType", "RangeFrame", "PyRanges")


class UnaryRangeFrameOperation[T: "RangeFrame"](Protocol):
    """A protocol for unary operations on RangeFrames."""

    def __call__(self, df: T, **kwargs: Any) -> "pd.DataFrame":
        """Perform the operation on the RangeFrame.

        Examples: cluster, merge, split, etc.
        """
        ...


class BinaryRangeFrameOperation[T: "RangeFrame"](Protocol):
    """A protocol for binary operations on RangeFrames."""

    def __call__(self, df: T, df2: T, **kwargs: Any) -> "pd.DataFrame":
        """Perform the operation on the pair of RangeFrames.

        Examples: overlap, nearest, join, etc.
        """
        ...

# Define the Literal type
VALID_OVERLAP_TYPE = Literal["first", "containment", "all"]

# Extract the options from the Literal type
VALID_OVERLAP_OPTIONS = list(get_args(VALID_OVERLAP_TYPE))
OVERLAP_FIRST, OVERLAP_CONTAINMENT, OVERLAP_ALL = VALID_OVERLAP_OPTIONS

CHROM_COL: Final = "Chromosome"
START_COL = "Start"
END_COL = "End"
STRAND_COL = "Strand"
RANGE_COLS = [START_COL, END_COL]
CHROM_AND_STRAND_COLS = [CHROM_COL, STRAND_COL]
GENOME_LOC_COLS = [CHROM_COL, *RANGE_COLS]
GENOME_LOC_COLS_WITH_STRAND = [*GENOME_LOC_COLS, STRAND_COL]

BIGWIG_SCORE_COL = "Score"
FRAME_COL = "Frame"

FORWARD_STRAND: Final = "+"
REVERSE_STRAND: Final = "-"
VALID_GENOMIC_STRAND_TYPE = Literal["+", "-"]
VALID_GENOMIC_STRAND_INFO = [FORWARD_STRAND, REVERSE_STRAND]

VALID_BY_OPTIONS = str | Iterable[str] | None

STRAND_AUTO: Final = "auto"
VALID_STRAND_TYPE = Literal["auto"] | bool
VALID_STRAND_OPTIONS = [STRAND_AUTO, True, False]

STRAND_BEHAVIOR_AUTO: Final = "auto"
STRAND_BEHAVIOR_SAME: Final = "same"
STRAND_BEHAVIOR_OPPOSITE: Final = "opposite"
STRAND_BEHAVIOR_IGNORE: Final = "ignore"
VALID_STRAND_BEHAVIOR_TYPE = Literal["auto", "same", "opposite", "ignore"]
VALID_STRAND_BEHAVIOR_OPTIONS = [
    STRAND_BEHAVIOR_SAME,
    STRAND_BEHAVIOR_OPPOSITE,
    STRAND_BEHAVIOR_AUTO,
    STRAND_BEHAVIOR_IGNORE,
]

JOIN_OUTER: Final = "outer"
JOIN_INNER: Final = "inner"
JOIN_RIGHT: Final = "right"
JOIN_LEFT: Final = "left"
VALID_JOIN_TYPE = Literal["inner", "left", "outer", "right"]
VALID_JOIN_OPTIONS = [JOIN_INNER, JOIN_LEFT, JOIN_OUTER, JOIN_RIGHT]

NEAREST_ANY_DIRECTION: Final = "any"
NEAREST_UPSTREAM: Final = "upstream"
NEAREST_DOWNSTREAM: Final = "downstream"
VALID_NEAREST_TYPE = Literal["any", "upstream", "downstream"]
VALID_NEAREST_OPTIONS = [
    NEAREST_ANY_DIRECTION,
    NEAREST_UPSTREAM,
    NEAREST_UPSTREAM,
    None,
]

TEMP_INDEX_COL = "__temp_index__"
TEMP_TYPE_COL = "__temp_type__"
TEMP_START_COL = "__temp_start__"
TEMP_END_COL = "__temp_end__"
TEMP_STRAND_COL = "__temp_strand__"
TEMP_NAME_COL = "__temp_name__"
TEMP_CUMSUM_COL = "__temp_cumsum__"
TEMP_NUM_COL = "__temp_num__"
TEMP_ID_COL = "__temp_id__"
TEMP_MIN_COL = "__temp_min__"
TEMP_MAX_COL = "__temp_max__"
TEMP_CLUSTER_COL = "__temp_cluster__"
TEMP_LENGTH_COL = "__temp_length__"
TEMP_START_SLACK_COL = "__temp_start_slack__"
TEMP_END_SLACK_COL = "__temp_end_slack__"

DEFAULT_CLUSTER_ID_COL = "Cluster"

VALID_BY_TYPES = str | Iterable[str] | None

PANDAS_COMPRESSION_TYPE = Literal["infer", "gzip", "bz2", "zip", "xz", "zstd"] | dict[str, Any] | None
