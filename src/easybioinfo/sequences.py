import itertools
__all__=['translate', 'reverse_complement']

##### reverse complement
# thanks to Devon Ryan at https://bioinformatics.stackexchange.com/questions/3583/what-is-the-fastest-way-to-get-the-reverse-complement-of-a-dna-sequence-in-pytho
complement =    str.maketrans("ACTGactg", "TGACtgac")
rnacomplement = str.maketrans("ACUGacug", "UGACugac")
def reverse_complement(seq, is_RNA=False):
    """ Reverse complement a DNA sequence

    Parameters
    ----------
    seq : str
        nucleotide sequence in DNA format (characters: ATGC)

    is_RNA : bool
        use this to provide the input seq in RNA format instead (characters: AUGC)

    Returns
    -------
    revcompseq : str
        reverse complement nucleotide sequence in DNA format (or RNA if is_RNA was set to True)

    Note
    ---- 
    Characters that are not upper or lowercase ATGC (or AUGC if is_RNA) are left unchanged
    """
    if not is_RNA:
        return seq.translate(complement)[::-1]
    else:
        return seq.translate(rnacomplement)[::-1]


#### Translate function
# build alternative genetic code translation tables based on NCBI codes
genetic_codes={} # will store dict of dicts: {genetic_code_id: {codon:aminoacid} }
genetic_codes_AAs={"1":"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "2":"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
                   "3":"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "4":"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "5":"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
                   "6":"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "9":"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                   "10":"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "11":"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "12":"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "13":"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
                   "14":"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                   "16":"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "21":"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                   "22":"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "23":"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "24":"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
                   "25":"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "26":"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "27":"FFLLSSSSYYQQCCWWLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "28":"FFLLSSSSYYQQCCWWLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "29":"FFLLSSSSYYYYCC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "30":"FFLLSSSSYYEECC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "31":"FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                   "33":"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
}

for gc_code in genetic_codes_AAs:
    genetic_codes[gc_code]={'---':'-'} 
    for codon_index, combo in enumerate( itertools.product( 'TCAG', repeat=3) ):   # iterates over ('T', 'T, 'T),  ('T', 'T', 'C'), ... etc 
        codon=''.join(combo)
        genetic_codes[gc_code][codon]=genetic_codes_AAs[gc_code][codon_index]
    genetic_codes[gc_code+'+U']=genetic_codes[gc_code].copy()
    genetic_codes[gc_code+'+U']['TGA']='U'
        
def translate(seq, genetic_code='1', unknown='X'):
    """ Translate a coding sequence into protein

    Parameters
    ----------
    seq : str
        nucleotide sequence in DNA format (characters: ATGC)
    genetic_code : str | dict
        string-converted NCBI index for genetic code (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
        or dictionary with keys for each codon, values are amino acids (remember to 
        include the translation for gaps ``'---':'-'``)
        or string-converted NCBI index with a '+U' suffix to have UGA as selenocysteine (U character)
    unknown : str | None
        codons that are not found in the genetic code table will be translated as this character
        if None, finding an unknown codon will raise an exception instead
    
    Returns
    -------
    pep : str
        protein sequence resulting from translation, with gaps as '-' and unknown characters as 'X'

    Warning
    -------
    This function expects uppercase DNA as input. 
    'U' or lowercase characters will result in 'X' characters as translation.
    To provide more flexible input, pre-process input 
    with ``seq.upper().replace('U', 'T')``

    """
    if genetic_code in genetic_codes:
        codon_table=genetic_codes[genetic_code]
    elif type(genetic_code) is dict:
        codon_table=genetic_code
    else:
        raise Exception(f'translate ERROR genetic_code input not recognized: {genetic_code}')

    output=[]
    for pos in range(0, len(seq), 3):
        codon=seq[pos:pos+3]
        if codon in codon_table:
            output.append( codon_table[codon] )
        elif unknown is not None:
            output.append( unknown )
        else:
            raise Exception(f'translate ERROR cannot find codon {codon} (pos {pos}-{pos+3}) in genetic_code!')
        
    return ''.join(output)


    
        
    
