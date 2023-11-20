from typing import Literal

CHROM_COL = "Chromosome"
START_COL = "Start"
END_COL = "End"
STRAND_COL = "Strand"
RANGE_COLS = [START_COL, END_COL]
GENOME_LOC_COLS = [CHROM_COL, *RANGE_COLS]
GENOME_LOC_COLS_WITH_STRAND = [*GENOME_LOC_COLS, STRAND_COL]

FORWARD_STRAND = "+"
REVERSE_STRAND = "-"
VALID_GENOMIC_STRAND_INFO = [FORWARD_STRAND, REVERSE_STRAND]

STRAND_AUTO = "auto"
VALID_STRAND_OPTIONS = Literal["auto"] | True | False

STRAND_BEHAVIOR_AUTO = "auto"
STRAND_BEHAVIOR_SAME = "same"
STRAND_BEHAVIOR_OPPOSITE = "opposite"
STRAND_BEHAVIOR_IGNORE = "ignore"
VALID_STRAND_BEHAVIOR_TYPE = (
    Literal["auto"] | Literal["same"] | Literal["opposite"] | Literal["ignore"]
)
VALID_STRAND_BEHAVIOR_OPTIONS = [
    STRAND_BEHAVIOR_SAME,
    STRAND_BEHAVIOR_OPPOSITE,
    STRAND_BEHAVIOR_AUTO,
    STRAND_BEHAVIOR_IGNORE,
]

OVERLAP_FIRST = "first"
OVERLAP_CONTAINMENT = "containment"
OVERLAP_ALL = "all"
VALID_OVERLAP_TYPE = Literal["first"] | Literal["containment"] | Literal["all"]
VALID_OVERLAP_OPTIONS = [OVERLAP_FIRST, OVERLAP_CONTAINMENT, OVERLAP_ALL]

JOIN_OUTER = "outer"
JOIN_INNER = "inner"
JOIN_RIGHT = "right"
JOIN_LEFT = "left"
VALID_JOIN_TYPE = (
    Literal["inner"] | Literal["left"] | Literal["outer"] | Literal["right"]
)
VALID_JOIN_OPTIONS = [JOIN_INNER, JOIN_LEFT, JOIN_OUTER, JOIN_RIGHT]

TEMP_INDEX_COL = "__temp_index__"
TEMP_START_COL = "__temp_start__"
TEMP_END_COL = "__temp_end__"
TEMP_NAME_COL = "__temp_name__"
