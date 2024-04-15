import itertools

import pyranges as pr  # noqa: F401

# thanks to Devon Ryan at https://bioinformatics.stackexchange.com/questions/3583/what-is-the-fastest-way-to-get-the-reverse-complement-of-a-dna-sequence-in-pytho
complement = str.maketrans("ACTGactg", "TGACtgac")
rnacomplement = str.maketrans("ACUGacug", "UGACugac")

# build alternative genetic code translation tables based on NCBI codes
GENETIC_CODES = {}  # will store dict of dicts: {genetic_code_id: {codon:aminoacid} }
GENETIC_CODE_AAS = {
    "1": "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "2": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    "3": "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "4": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "5": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    "6": "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "9": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "10": "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "11": "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "12": "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "13": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    "14": "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "16": "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "21": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "22": "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "23": "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "24": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
    "25": "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "26": "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "27": "FFLLSSSSYYQQCCWWLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "28": "FFLLSSSSYYQQCCWWLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "29": "FFLLSSSSYYYYCC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "30": "FFLLSSSSYYEECC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "31": "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "33": "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
}

for gc_code in GENETIC_CODE_AAS:
    GENETIC_CODES[gc_code] = {"---": "-"}
    for codon_index, combo in enumerate(
        itertools.product("TCAG", repeat=3),
    ):  # iterates over ('T', 'T, 'T),  ('T', 'T', 'C'), ... etc
        codon = "".join(combo)
        GENETIC_CODES[gc_code][codon] = GENETIC_CODE_AAS[gc_code][codon_index]
    GENETIC_CODES[gc_code + "+U"] = GENETIC_CODES[gc_code].copy()
    GENETIC_CODES[gc_code + "+U"]["TGA"] = "U"


def reverse_complement(seq: str, *, is_rna: bool = False, check: bool = False) -> str:
    """Reverse complement a DNA sequence.

    The input sequence is in DNA alphabet, upper or lowercase (ATGCatgc). The case is maintained in output.
    Other characters are left unchanged.

    Parameters
    ----------
    seq : str
        nucleotide sequence in DNA format (characters: ATGCatgc)

    is_rna : bool
        use this to provide the input seq in RNA format instead (characters: AUGCaugc)

    check : bool
        check if the input string contains only characters present in the translation table, raising an error if not

    Returns
    -------
    revcompseq : str
        reverse complement nucleotide sequence in DNA format (or RNA if is_rna was set to True)

    Examples
    --------
    >>> pr.seqs.reverse_complement("ATGAAATTTGGGTGA")
    'TCACCCAAATTTCAT'

    >>> pr.seqs.reverse_complement("AUGAAAUUUGGGUGA", is_rna=True)
    'UCACCCAAAUUUCAU'

    >>> pr.seqs.reverse_complement("aaaATCcccGGG")
    'CCCgggGATttt'

    >>> pr.seqs.reverse_complement("ATCWWWCCCTTT")
    'AAAGGGWWWGAT'

    >>> pr.seqs.reverse_complement("AUGAAATGGGTGA", check=True)
    Traceback (most recent call last):
    ...
    ValueError: One or more characters in the input string are not present in the translation table.

    """
    transtable = complement if not is_rna else rnacomplement
    if check and not all(char in transtable for char in set(seq)):
        msg = "One or more characters in the input string are not present in the translation table."
        raise ValueError(msg)
    return seq.translate(transtable)[::-1] if not is_rna else seq.translate(transtable)[::-1]


GENETIC_CODE_OPTION = str | int | dict


def translate(
    seq: str,
    *,
    genetic_code: GENETIC_CODE_OPTION = "1",
    unknown: str = "X",
    sanitize: bool = True,
    cache: int | bool = False,
) -> str:
    """Translate a coding sequence into protein.

    The input sequence is expected to be a DNA sequence in uppercase format (characters: ATGC).
    Incomplete codons at the end of the sequence, as well as non-canonical codons, result in the unknown character "X".

    Parameters
    ----------
    seq : str
        nucleotide sequence in uppercase DNA format (characters: ATGC)

    genetic_code : int | str | dict
        int or string-converted NCBI index for genetic code (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
        or dictionary with keys for each codon, values are amino acids (remember to
        include the translation for gaps ``'---':'-'``)
        or string-converted NCBI index with a '+U' suffix to have UGA as selenocysteine (U character)

    unknown : str | None
        codons that are not found in the genetic code table will be translated as this character
        if None, finding an unknown codon will raise an exception instead

    sanitize : bool
        Whether the input is converted to DNA uppercase (from lowercase and/or RNA) to ensure translation works.
        Set to False if you are sure sequences are already DNA uppercase, to speed up execution

    cache : bool | int
        speeds up translation by caching the result of translation of multicodon strings.
        The 1st time, the function is slow since precomputing all results; then, it is ~3 times faster than non-caching translate.
        With cache=True, all 3-codon translations are cached (memory 10Mb, precompute ~ 230ms).
        Provide an int to define how many codons to cache; this is approx the speedup that will be obtained.
        Note: memory and precomputing grow exponentially with the N of codons cached;
        You may use pyranges.seqs.clear_kmer_memory() to free memory.

    Returns
    -------
    pep : str
        protein sequence resulting from translation, with gaps as '-' and unknown characters as 'X'

    Examples
    --------
    >>> pr.seqs.translate("ATGAAATTTGGGTGA")
    'MKFG*'

    >>> pr.seqs.translate("ATGTTGCTGAA")
    'MLLX'

    # translate with vertebrate mithochondrial genetic code
    >>> pr.seqs.translate("ATGAAATTTGGGTGA", genetic_code=2)
    'MKFGW'

    # translate with a custom genetic code (all codons starting with A are translated as A, the rest as Q)
    >>> gc={codon:'A' if codon.startswith('A') else 'Q' for codon in map(''.join, itertools.product("TCAG", repeat=3))}
    >>> pr.seqs.translate("ATGAAATTTGGGTGA", genetic_code=gc)
    'AAQQQ'

    >>> pr.seqs.translate("AUGAAATTtGGGTGA", sanitize=False)
    'XKXG*'

    >>> pr.seqs.translate("AUGAAATTtGGGTGA", sanitize=False, unknown=None)
    Traceback (most recent call last):
    ...
    ValueError: translate ERROR cannot find codon AUG (pos 0-3) in genetic_code!


    # caching makes sense when translating many sequences
    # for large enough data, translate with cache=True is 3x faster than non-caching translate
    >>> import random
    >>> random.seed(42)
    >>> many_seqs = [ "".join(random.choices("ATGC", k=1000)) for _ in range(100000) ]
    >>> translated_seqs = [pr.seqs.translate(seq, cache=True) for seq in many_seqs]
    >>> translated_seqs[0][:50]
    'DRHKKEHVLTRLRGINVIGDRANYLS*CYESYRQGDS*PGRDTLPYHV*D'

    """
    if isinstance(genetic_code, dict):
        codon_table = genetic_code
    elif str(genetic_code) in GENETIC_CODES:
        codon_table = GENETIC_CODES[str(genetic_code)]
    else:
        msg = f"translate ERROR genetic_code input not recognized: {genetic_code}"
        raise ValueError(msg)

    if sanitize:
        seq = seq.upper().replace("U", "T")

    output = []
    if not cache:
        for pos in range(0, len(seq), 3):
            codon = seq[pos : pos + 3]
            if codon in codon_table:
                output.append(codon_table[codon])
            elif unknown is not None:
                output.append(unknown)
            else:
                msg = f"translate ERROR cannot find codon {codon} (pos {pos}-{pos+3}) in genetic_code!"
                raise ValueError(msg)
    else:
        cache = 3 if cache is True else cache
        kmer_codon_table = _get_kmer_codon_table(genetic_code, cache)
        step = cache * 3
        for pos in range(0, len(seq), step):
            multicodon = seq[pos : pos + step]
            output.append(
                kmer_codon_table[multicodon]
                if multicodon in kmer_codon_table
                # for final bits of sequences which are shorter than k codons, and also for sequences containing Ns:
                else translate(multicodon, genetic_code=genetic_code, unknown=unknown, cache=False),
            )

    return "".join(output)


kmer_codon_tables = {}


def clear_kmer_memory() -> None:
    """Clear memory used by the translation cache."""
    kmer_codon_tables.clear()


def _get_kmer_codon_table(genetic_code, k) -> dict[str, str]:
    if (genetic_code, k) not in kmer_codon_tables:
        all_codons = list(GENETIC_CODES["1"].keys())
        codon_table_multik = {
            multicodon: translate(multicodon, genetic_code=genetic_code)
            for multicodon in map("".join, itertools.product(*(all_codons for i in range(k))))
        }
        kmer_codon_tables[(genetic_code, k)] = codon_table_multik

    return kmer_codon_tables[(genetic_code, k)]
