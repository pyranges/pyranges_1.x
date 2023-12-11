from typing import TYPE_CHECKING

from pyranges.names import END_COL, REVERSE_STRAND, START_COL, TEMP_CUMSUM_COL, TEMP_INDEX_COL, TEMP_LENGTH_COL

if TYPE_CHECKING:
    import pyranges as pr


TOTAL_EXON_LENGTH_COL = "__temp_total_exon_len__"


def _spliced_subseq(
    scdf: "pr.PyRanges",
    start: int = 0,
    end: int | None = None,
    **kwargs,
) -> "pr.PyRanges":
    original_class = scdf.__class__
    if scdf.empty:
        return scdf

    scdf = scdf.copy()
    orig_order = scdf.index.copy()

    by_argument_given = kwargs.get("by")
    _by = kwargs.get("by", TEMP_INDEX_COL)
    by = [_by] if isinstance(_by, str) else (_by or [])

    strand = kwargs.get("strand")

    # at this point, strand is False if 1. spliced_subsequence was called with strand=False or
    #                                   2. it was called with strand=None and self is not stranded
    # or it can be '+' or '-' if:
    #  1. it was input as True to spliced_subsequence and passed  to pyrange_apply_single as True,
    #     which updates it to '-' or '+' before calling _spliced_subseq, or
    #  2. it was called with strand=None and self is stranded

    if strand and not scdf.strand_values_valid:
        msg = "Cannot have strand=True on unstranded pyranges!"
        raise AssertionError(msg)

    scdf.insert(scdf.shape[1], TEMP_LENGTH_COL, scdf.End - scdf.Start)
    scdf.insert(scdf.shape[1], TEMP_INDEX_COL, scdf.index)

    g = scdf.groupby(by, dropna=False)
    scdf.insert(scdf.shape[1], TEMP_CUMSUM_COL, g[TEMP_LENGTH_COL].cumsum())

    end = scdf[TEMP_CUMSUM_COL].max() if end is None else end

    minstart_idx = g[TEMP_INDEX_COL].first()

    if start < 0 or (end is not None and end < 0):
        # len_per_transc is total sum of exon length per transcript
        len_per_transc = scdf.loc[g[TEMP_INDEX_COL].last(), [*by, TEMP_CUMSUM_COL]].rename(
            columns={TEMP_CUMSUM_COL: TOTAL_EXON_LENGTH_COL},
        )

        # exp_len_per_transc has same rows of scdf with total sum of exon length
        # had to add bits to keep the order of rows right, or merge would destroy it
        if by_argument_given:
            exp_len_per_transc = (
                scdf.loc[:, [*by, TEMP_INDEX_COL]]
                .merge(len_per_transc, on=by)
                .set_index(TEMP_INDEX_COL)
                .loc[scdf.index]
            )
        else:
            exp_len_per_transc = scdf.loc[:, by].merge(len_per_transc, on=by).set_index(TEMP_INDEX_COL).loc[scdf.index]

        if start < 0:
            start = exp_len_per_transc[TOTAL_EXON_LENGTH_COL] + start

        if end is not None and end < 0:
            end = exp_len_per_transc[TOTAL_EXON_LENGTH_COL] + end

    cs_start = g[TEMP_CUMSUM_COL].shift(1, fill_value=0)
    cs_start.loc[minstart_idx] = 0

    cs_end = scdf[TEMP_CUMSUM_COL]

    # NOTE
    # here below, start is a scalar if originally provided > 0, or a Series if < 0
    #             end is a scalar if originally None or provided >0, or a Series if provided < 0
    if strand == REVERSE_STRAND:  # and use_strand:
        start_adjustments = start - cs_start
        adjust_start = start_adjustments > 0
        scdf.loc[adjust_start, END_COL] -= start_adjustments[adjust_start].astype(scdf[END_COL].dtype)

        end_adjustments = cs_end - end
        adjust_end = end_adjustments > 0
        scdf.loc[adjust_end, START_COL] += end_adjustments[adjust_end].astype(scdf[START_COL].dtype)
    else:
        start_adjustments = start - cs_start
        adjust_start = start_adjustments > 0
        scdf.loc[adjust_start, START_COL] += start_adjustments[adjust_start].astype(scdf[START_COL].dtype)

        end_adjustments = cs_end - end
        adjust_end = end_adjustments > 0
        scdf.loc[adjust_end, END_COL] -= end_adjustments[adjust_end].astype(scdf[END_COL].dtype)

    _scdf = scdf.loc[orig_order]
    _scdf = _scdf[(_scdf[START_COL] < _scdf[END_COL])]

    return original_class(_scdf.drop([TEMP_INDEX_COL, TEMP_LENGTH_COL, TEMP_CUMSUM_COL], axis=1))
