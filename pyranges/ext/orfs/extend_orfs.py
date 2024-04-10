import logging
import sys
import warnings
from collections.abc import Iterable
from typing import Any, Literal

import numpy as np
import pandas as pd

import pyranges as pr
from pyranges.core.names import CHROM_COL, END_COL, FORWARD_STRAND, REVERSE_STRAND, START_COL, STRAND_COL
from pyranges.core.pyranges_helpers import mypy_ensure_pyranges

## Joan Pallares started this function
#  Marco Mariotti later optimized it, made it independent from pyranges

logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

DIRECTION_OPTIONS = Literal["up", "down"]

STARTS_NUCLEOTIDE_SEQ = ["ATG"]
STOPS_NUCLEOTIDE_SEQ = ["TAG", "TGA", "TAA"]

N_PRINTED_NON_MULTIPLE_3 = 10


def extend_orfs(  # noqa: C901,PLR0912,PLR0915
    self: pr.PyRanges,
    fasta_path: str,
    transcript_id: str | list[str] | None = None,
    *,
    starts: list[str] = STARTS_NUCLEOTIDE_SEQ,
    stops: list[str] = STOPS_NUCLEOTIDE_SEQ,
    keep_off_bounds: bool = False,
    direction: Iterable[DIRECTION_OPTIONS] | DIRECTION_OPTIONS | None = None,
    chunk_size: int = 900,
    record_extensions: bool = False,
    verbose: bool = False,
) -> pr.PyRanges:
    """Extend PyRanges intervals to form complete open reading frames.

    The input intervals are extended their next Stop codon downstream, and to their leftmost Start codon upstream
    before encountering a Stop.

    Parameters
    ----------
    self : PyRanges
        Input CDS intervals.

    fasta_path : location of the Fasta file from which the sequences
        for the extensions will be retrieved.

    transcript_id : str or list of str or None
        Name(s) of column(s) to group intervals into transcripts

    starts : list containing the nucleotide pattern to look for upstream.
        if not provided, ORFs are delimited by stops
        Default ['ATG']

    stops : list containing the nucleotide pattern to look for downstream.
        Default ['TAG', 'TGA', 'TAA']

    keep_off_bounds : if True, those intervals that reached out of bounds during extension
        without finding any stop are returned in their largest (3-nt multiple) extension.
        In this case, these intervals will not begin with a start or end with a stop

    direction : whether the extension should be upstream ('up'), downstream
        ('down') or both. Default (None) means: ['up', 'down']

    chunk_size : the amount of nucleotides to be extended on each iteration.
        Default 900.

    record_extensions : if True, add columns extension_up and extension_down
        with the extensions amounts. Default: False

    verbose : if True, print information about the progress of the extension.
        Default: False

    """
    try:
        import pyfaidx  # type: ignore[import]
    except ImportError:
        LOGGER.exception(
            "To use extend_orfs, pyfaidx must be installed. Use `conda install -c bioconda pyfaidx` or `pip install pyfaidx` to install pyfaidx.",
        )
        sys.exit(1)

    def pverbose(msg: Any) -> None:
        if verbose:
            LOGGER.info(msg)

    # Sanity Checks
    if (chunk_size % 3) != 0:
        msg = "Chunk size must be a multiple of three."
        raise AssertionError(msg)

    if not all(len(pattern) == len(STARTS_NUCLEOTIDE_SEQ[0]) for pattern in starts + stops):
        msg = "Ensure that all patterns have a length of 3 nt."
        raise AssertionError(msg)

    direction = direction if direction is not None else ["up", "down"]

    p = self.copy()

    if not p.strand_valid:
        msg = "Intervals must be have valid strands to call extend_orfs"
        raise AssertionError(msg)

    if transcript_id is None:
        cds_id = "Custom_ID"
        p[cds_id] = np.arange(len(p))  # Generates a New Column
    else:
        cds_id = transcript_id

    if isinstance(cds_id, str):  # ensure it is a list
        cds_id = [cds_id]

    ##################
    # Show a warning if some transcript lengths are not divisible by 3:
    p["__length"] = p[END_COL] - p[START_COL]
    module_3 = p.groupby(by=cds_id).__length.sum().mod(3)  # noqa: SLF001
    ndv = module_3[module_3 != 0].index
    if len(ndv) > 0:
        warnings.warn("\nWARNING! Some input lenghts are not divisible by 3:", stacklevel=2)
        warnings.warn(
            " ".join(ndv[:N_PRINTED_NON_MULTIPLE_3]) + (" and others " if len(ndv) > N_PRINTED_NON_MULTIPLE_3 else ""),
            stacklevel=2,
        )

    p = p.drop_and_return(["__length"], axis=1)
    ##################

    # Load Sequence Data from a Fasta file
    fs = pyfaidx.Fasta(fasta_path)  # pyfaidx_fasta object, fasta sequences

    # get a minimal interval per group, with min start and max end
    minp = (
        (p[[STRAND_COL, START_COL, END_COL, CHROM_COL]+cds_id])
        .groupby(cds_id, as_index=False)
        .agg(
            {
                STRAND_COL: "first",
                START_COL: "min",
                END_COL: "max",
                CHROM_COL: "first",
            },
        )
    )
    ncds = len(minp)
    nexons = len(p)

    # add chromosome sizes as column
    minp = minp.merge(
        pd.DataFrame(
            {
                CHROM_COL: [chrom for chrom, _ in fs.items()],
                "__chromsize": [len(fs[chrom]) for chrom, _ in fs.items()],
            },
        ),
        on=CHROM_COL,
    ).set_index(cds_id)

    if len(minp) != ncds:
        msg = "ERROR some sequences were not found in the pyfaidx object"
        raise AssertionError(msg)

    ######################  Extend Sequence Upstream ###########################

    if "up" in direction or direction == "up":
        pup = minp.copy()
        pup["__out_of_bounds"] = False
        pup["__up_stop_dist"] = -1  # -1 means not found
        pup["__up_start_dist"] = -1  # -1 means not found
        selector = pd.Series(
            data=True,
            index=pup.index,
        )  # initialized as all True; keeps track of intervals with work to do
        ic = 0  # iteration counter

        while (selector).any():
            pverbose(
                f"Upstream | iteration {ic},  left to do: {selector.sum()}",
            )
            pverbose(pup[selector])
            # selector= ... everything which is not done nor out of bounds
            _extend_df(pup, "5", chunk_size, selector)
            _subseq_df(pup, "5", chunk_size, selector)
            _correct_bounds(pup, adjust="5")
            _get_seqs(pup, fs, selector)  # seqs length is only selector

            ## looking for a stop if didn't find it already
            z = pup.loc[selector, "__seq"].apply(lambda x: _find_rightmost_stop(x, stops))
            found_stop_selector = selector.copy()  # made to index pup
            found_stop_selector[found_stop_selector] = z != -1
            pup.loc[found_stop_selector, "__up_stop_dist"] = z[z != -1] + ic * chunk_size

            if keep_off_bounds:
                # focus on intervals that reached out of bounds in this iteration and did not find a stop:
                #  record their max extension as up_stop_dist and __up_start_dist so they are considered whether starts was provided or not
                # btw this means these intervals can be identified by looking at rows with __up_start_dist == __up_stop_dist
                no_stop_off_bounds = found_stop_selector.copy()
                no_stop_off_bounds[no_stop_off_bounds] = z == -1
                pup.loc[no_stop_off_bounds, "__up_stop_dist"] = (
                    pup.loc[no_stop_off_bounds, END_COL] - pup.loc[no_stop_off_bounds, START_COL]
                ) + ic * chunk_size
                pup.loc[no_stop_off_bounds, "__up_start_dist"] = pup.loc[no_stop_off_bounds, "__up_stop_dist"]

            ## looking for a start
            if starts:
                z = pup.loc[selector, ["__seq", "__up_stop_dist"]].apply(
                    lambda x: _find_leftmost_start(x, starts, chunk_size),
                    axis=1,
                )
                ## if start is not found now, but was found in a previous iteration, the previous will be kept (not overwriting a -1)
                found_start_selector = selector.copy()  # made to index pup
                found_start_selector[found_start_selector] = z != -1
                pup.loc[found_start_selector, "__up_start_dist"] = z[z != -1] + ic * chunk_size

            pup.loc[found_stop_selector, "__seq"] = ""  ## frees memory
            pverbose(" ---> After iteration:")
            pverbose(pup[selector])

            selector = (pup.__up_stop_dist == -1) & ~(pup.__out_of_bounds)  # noqa: SLF001
            ic += 1

        ext_up = (
            pup["__up_stop_dist"].rename("extension_up")
            if not starts
            else pup["__up_start_dist"].rename("extension_up")
        )
        ext_up[ext_up == -1] = 0
        del pup

    else:  # No extension upstream (i.e. the extension is 0)
        ext_up = pd.Series(0, index=minp.index)

    ######################  Extend Sequence Downstream #########################
    if "down" in direction or direction == "down":
        pdo = minp.copy()
        pdo["__out_of_bounds"] = False
        pdo["__up_stop_dist"] = -1  # -1 means not found
        selector = pd.Series(
            data=True,
            index=pdo.index,
        )  # initialized as all True; keeps track of intervals with work to do
        ic = 0  # iteration counter

        while (selector).any():
            pverbose(
                f"Downstream | iteration {ic},  left to do: {selector.sum()}",
            )
            pverbose(pdo[selector])
            # selector= ... everything which is not done nor out of bounds
            _extend_df(pdo, "3", chunk_size, selector)
            _subseq_df(pdo, "3", chunk_size, selector)
            _correct_bounds(pdo, adjust="3")
            _get_seqs(pdo, fs, selector)  # seqs length is only selector

            ## looking for a stop if didn't find it already
            z = pdo.loc[selector, "__seq"].apply(lambda x: _find_leftmost_stop(x, stops))
            found_stop_selector = selector.copy()  # made to index pdo
            found_stop_selector[found_stop_selector] = z != -1
            pdo.loc[found_stop_selector, "__up_stop_dist"] = z[z != -1] + ic * chunk_size

            if keep_off_bounds:
                # focus on intervals that reached out of bounds in this iteration and did not find a stop:
                #  record their max extension as up_stop_dist
                no_stop_off_bounds = found_stop_selector.copy()
                no_stop_off_bounds[no_stop_off_bounds] = z == -1
                pdo.loc[no_stop_off_bounds, "__up_stop_dist"] = (
                    pdo.loc[no_stop_off_bounds, END_COL] - pdo.loc[no_stop_off_bounds, START_COL]
                ) + ic * chunk_size

            pdo.loc[found_stop_selector, "__seq"] = ""  ## frees memory
            pverbose(" ---> After iteration:")
            pverbose(pdo[selector])
            selector = (pdo.__up_stop_dist == -1) & ~(pdo.__out_of_bounds)  # noqa: SLF001
            ic += 1

        ext_down = pdo["__up_stop_dist"].rename("extension_down")
        ext_down[ext_down == -1] = 0
        del pdo

    else:  # No extension downstream (i.e. the extension is 0)
        ext_down = pd.Series(0, index=minp.index)

    ### Extensions have been determined. Now let's apply them to the original DF
    p["__order"] = np.arange(len(p))
    zp = p.merge(
        pd.DataFrame({"__extension_up": ext_up, "__extension_down": ext_down}),
        left_on=cds_id,
        right_index=True,
    )
    if len(zp) != nexons:
        msg = "ERROR malformed gene structures"
        raise AssertionError(msg)

    p = p.drop(columns=["__order"])
    p = zp

    p=p.sort_ranges()
    _extend_groups(p, cds_id)
    p = mypy_ensure_pyranges(p.sort_values("__order"))  # restore order

    # getting ready to return
    if record_extensions:
        p = mypy_ensure_pyranges(
            p.rename(
                columns={
                    "__extension_up": "extension_up",
                    "__extension_down": "extension_down",
                },
            ),
        )
        to_drop = ["__order"]
    else:
        to_drop = ["__extension_up", "__extension_down", "__order"]

    if transcript_id is None:
        to_drop.extend(cds_id)

    return mypy_ensure_pyranges(p.drop_and_return(to_drop, axis=1))


# MM # efficients methods to work inplace. All methods with df as input have one row per CDS group
def _extend_df(df, direction, extension, selector) -> None:
    """Inplace."""
    if direction == "5":
        df.loc[selector & (df.Strand == FORWARD_STRAND), [START_COL]] -= extension
        df.loc[selector & (df.Strand == "-"), [END_COL]] += extension
    elif direction == "3":
        df.loc[selector & (df.Strand == REVERSE_STRAND), [START_COL]] -= extension
        df.loc[selector & (df.Strand == FORWARD_STRAND), [END_COL]] += extension


def _extend_groups(dfg, cds_id) -> None:
    """In place, dfg may contain multiple exons."""
    first_exon_indices = dfg.groupby(cds_id).apply(lambda x: x.index[0]).reset_index(drop=True)
    last_exon_indices = dfg.groupby(cds_id).apply(lambda x: x.index[-1]).reset_index(drop=True)

    # Convert first_exon_indices to a set for faster intersection operations
    first_exon_index_set = set(first_exon_indices)
    last_exon_index_set = set(last_exon_indices)

    sp = dfg.Strand == FORWARD_STRAND
    sm = dfg.Strand == REVERSE_STRAND
    # Get indices where sp is True
    positive_strand_indices = set(sp.index[sp])
    negative_strand_indices = set(sm.index[sm])

    # Find the intersection of first exon indices and positive strand indices
    first_exon_positive_strand_indices = list(first_exon_index_set & positive_strand_indices)
    last_exon_positive_strand_indices = list(last_exon_index_set & positive_strand_indices)

    # same for negative strand
    first_exon_negative_strand_indices = list(first_exon_index_set & negative_strand_indices)
    last_exon_negative_strand_indices = list(last_exon_index_set & negative_strand_indices)

    #print(dfg, first_exon_indices, sp)

    dfg.loc[first_exon_positive_strand_indices, START_COL] -= dfg.loc[first_exon_positive_strand_indices, "__extension_up"]
    dfg.loc[last_exon_positive_strand_indices, END_COL] += dfg.loc[last_exon_positive_strand_indices, "__extension_down"]
    dfg.loc[first_exon_negative_strand_indices, END_COL] += dfg.loc[first_exon_negative_strand_indices, "__extension_up"]
    dfg.loc[last_exon_negative_strand_indices, START_COL] -= dfg.loc[last_exon_negative_strand_indices, "__extension_down"]


    # dfg.loc[sp[first_exon_indices], START_COL] -= dfg.loc[sp[first_exon_indices], "__extension_up"]
    # dfg.loc[sp[last_exon_indices], END_COL] += dfg.loc[sp[last_exon_indices], "__extension_down"]
    # dfg.loc[sm[first_exon_indices], END_COL] += dfg.loc[sm[first_exon_indices], "__extension_up"]
    # dfg.loc[sm[last_exon_indices], START_COL] -= dfg.loc[sm[last_exon_indices], "__extension_down"]


def _get_seqs(df, fs, selector) -> None:
    """Creates/updates column __seq."""
    seqs = []
    d = df[selector]
    for start, end, strand, chrom in zip(d.Start, d.End, d.Strand, d.Chromosome, strict=False):
        _fasta = fs[chrom]
        if strand == REVERSE_STRAND:
            seqs.append((-_fasta[start:end]).seq.upper())  # reverse complement
        elif strand == FORWARD_STRAND:
            seqs.append(_fasta[start:end].seq.upper())
    df.loc[selector, "__seq"] = seqs


def _correct_bounds(df, adjust=None) -> None:
    """Similar to pyranges genome_bounds, inplace, but see below.

    but also: adjust the 5' or 3' end of intervals so they are multiples of 3
    and set boolean attribute __out_of_bounds.
    """
    out_end = df[END_COL] > df["__chromsize"]
    df.loc[out_end, END_COL] = df.loc[out_end, [END_COL, "__chromsize"]].min(axis=1)
    out_start = df[START_COL] < 0
    df.loc[out_start, START_COL] = 0

    sp = df.Strand == FORWARD_STRAND
    sm = df.Strand == REVERSE_STRAND

    if adjust == "5":
        df.loc[out_start & sp, START_COL] += (df.loc[out_start & sp, END_COL] - df.loc[out_start & sp, START_COL]).mod(
            3,
        )
        df.loc[out_end & sm, END_COL] -= (df.loc[out_end & sm, END_COL] - df.loc[out_end & sm, START_COL]).mod(3)

    elif adjust == "3":
        df.loc[out_end & sp, END_COL] -= (df.loc[out_end & sp, END_COL] - df.loc[out_end & sp, START_COL]).mod(3)
        df.loc[out_start & sm, START_COL] += (df.loc[out_start & sm, END_COL] - df.loc[out_start & sm, START_COL]).mod(
            3,
        )

    df.loc[(out_end | out_start), "__out_of_bounds"] = True


def _subseq_df(df, from_end, amount, selector) -> None:
    """Similar to pyranges subseq but inplace."""
    sp = selector & (df.Strand == FORWARD_STRAND)
    sm = selector & (df.Strand == REVERSE_STRAND)
    if from_end == "5":
        df.loc[sp, END_COL] = df.loc[sp, START_COL] + amount
        df.loc[sm, START_COL] = df.loc[sm, END_COL] - amount
    elif from_end == "3":
        df.loc[sp, START_COL] = df.loc[sp, END_COL] - amount
        df.loc[sm, END_COL] = df.loc[sm, START_COL] + amount


def _find_rightmost_stop(seq, motifs) -> int:
    """Return the distance from end of seq of the 1st codon position of the rightmost stop; e.g. 3 is for last codon."""
    for i in range(len(seq) - 3, -1, -3):
        if seq[i : i + 3] in motifs:
            return len(seq) - i
    return -1


def _find_leftmost_start(row, motifs, chunk_size) -> int:
    """Return distance from end of seq of 1st codon position of leftmost start AFTER the pos in column __up_stop_dist."""
    seq = row["__seq"]
    stop_pos = row["__up_stop_dist"]
    start_search = (3 + len(seq) - (stop_pos % chunk_size)) if stop_pos != -1 else 0
    # we search starting from the position next to the stop
    for i in range(start_search, len(seq), 3):
        if seq[i : i + 3] in motifs:
            return len(seq) - i
    return -1


def _find_leftmost_stop(seq, motifs) -> int:
    """Return the position (distance from the start of the string) of the 3rd codon position of the leftmost stop.

    e.g. 3 is for the first codon.
    """
    for i in range(0, len(seq), 3):
        if seq[i : i + 3] in motifs:
            return i + 3
    return -1

