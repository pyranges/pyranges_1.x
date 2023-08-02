import itertools
import warnings

__all__=['translate', 'reverse_complement', 'count_coding_sites', 'count_coding_changes', 'count_unique_changes']

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
   
        
    
codon2sitecount={}  # store computed data per codon
def count_coding_sites(cds, silent=False, cpg=False, nonsense=False):
    """ Counts the number of synonymous and non-synonymous sites for a input nucletoide coding sequence. 
  
    This function is instrumental for any dN/dS analysis. 
    As a single site can be partly non-syn and partly syn, the numbers returned are float (always multiple of one third).
    Gaps "-", if present, are removed silently.
    Codons with any character different from ACTG (for example, N) are skipped and a message is printed to stderr.
    
    Parameters
    ----------
    cds : str
        coding sequence in DNA format (characters: ATGC)
  
    silent : bool
        do not print warning messages for non-ATGC characters

    cpg : bool
        if True, separate counts for CpG sites only are provided, altering the structure of the returned tuple (see below).
        Note that the total counts are still provided. To obtain the number of nonCpG changes, subtract the number of CpG sites from the total number

    nonsense : bool
        if True, non-sense mutations (i.e. stop codons or mutations to a stop codon) are not counted as non-synonymous, and are returned separately, 
        altering the structure of the returned tuple (see below).

    Returns
    -------
    counts : tuple of floats

        * With default options, counts contains:                  (nonsyn, syn)
        * If cpg is True (but nonsense is not), counts contains:  (nonsyn, syn, cpg_nonsyn, cpg_syn)
        * If nonsense is True (but cpg is not), counts contains:  (nonsyn, syn, nonsense)
        * If both cpg and nonsense are true, counts contains:     (nonsyn, syn, nonsense, cpg_nonsyn, cpg_syn, cpg_nonsense)

    """   
    global codon2sitecount
    cds=cds.replace('-', '').upper()
    # these will result to be three times as much the actual values: I divide them as the very last step!
    syn=0
    nonsyn=0   
    cpg_syn=0
    cpg_nonsyn=0
    nonsensen=0
    cpg_nonsense=0
    dna_letts=set('ACTG')

    if len(cds)%3!=0:
        raise Exception("count_coding_sites ERROR the sequence must be composed of codons (length multiple of 3)")

    for i_codon in range(len(cds)//3):
      ## cycling codons
        codon=cds[i_codon*3:i_codon*3+3]
        if not all([lett in dna_letts  for lett in codon] ):       
            if not silent:
                warnings.warn(f'count_coding_sites WARNING skipping codon n.{i_codon+1} : {codon}')
            continue

        if codon in codon2sitecount and not cpg:
            # retrieving precomputed counts for this codon. Possible only if cpg don't need to be split
            n,s,x =codon2sitecount[codon]
        else:
            n,s,x=0,0,0
            if cpg:
                cn,cs,cx=0,0,0

            for i_within_codon in range(3):
                nt= codon[i_within_codon]
                i_cds=i_codon*3+i_within_codon
                syn_this_pos=0
                nonsense_this_pos=0
                for alt_nt in 'ACTG':
                    if alt_nt==nt:
                        continue
                    alt_codon=   codon[:i_within_codon]+alt_nt+codon[i_within_codon+1:]
                    if translate(alt_codon)==translate(codon):          
                        syn_this_pos+=1
                    elif "*" in (translate(alt_codon), translate(codon)):
                        nonsense_this_pos+=1
                nonsyn_this_pos=3-syn_this_pos-nonsense_this_pos

                if cpg: 
                    is_cpg= (   nt == 'G' and (  (  i_cds+1<len(cds)  and   cds[i_cds+1] =='C' )  or (i_cds!=0 and cds[i_cds-1] =='C' )    )     ) or \
                      (   nt == 'C' and (  (  i_cds+1<len(cds)  and   cds[i_cds+1] =='G' )  or (i_cds!=0 and cds[i_cds-1] =='G' )    )     )        # # G and the next or previous is C OR #C and the
                    if is_cpg:    
                        cs+= syn_this_pos
                        cn+= nonsyn_this_pos
                        cx+= nonsense_this_pos

                s+=syn_this_pos      
                n+=nonsyn_this_pos
                x+=nonsense_this_pos

            codon2sitecount[codon]=(n,s,x)

        syn+=s
        nonsyn+=n
        nonsensen+=x
        if cpg:
            cpg_syn+=cs
            cpg_nonsyn+=cn
            cpg_nonsense+=cx

    if cpg:
        if nonsense:
            return ( nonsyn/3, syn/3, nonsensen/3, cpg_nonsyn/3,  cpg_syn/3, cpg_nonsense/3 )
        else:
            return ( (nonsyn+nonsensen)/3, syn/3, (cpg_nonsyn+cpg_nonsense)/3,  cpg_syn/3 )
    elif nonsense:
        return ( nonsyn/3, syn/3, nonsensen/3 )
    else:
        ## default return value
        return ( (nonsyn+nonsensen)/3, syn/3 )

    

def count_coding_changes(cds, cds2, silent=True, cpg=False, nonsense=False):
    """Counts the number of synonymous and non-synonymous changes between two sequences.

       This function is instrumental for any dN/dS analysis.
       Sequences should be aligned using gaps "-". Gapped positions are not counted.
       Codons with any character different from ACTG (for example, N) are skipped and a message is printed to stderr.

      silent : bool, defaults to True
          do not print warning messages for non-ATGC characters

      cpg : bool
          if True, separate counts for CpG sites only are provided, altering the structure of the returned tuple (see below).
          Note that only cds is inspected to identity CpG sites, so with cpg=True the output may be different for (cds, cds2) and (cds2, cds)
          Note that the total counts are still provided. To obtain the number of nonCpG changes, subtract the number of CpG sites from the total number

      nonsense : bool
          if True, non-sense mutations (i.e. stop codons or mutations to a stop codon) are not counted as non-synonymous, and are returned separately, 
          altering the structure of the returned tuple (see below).

      Returns
      -------
      counts : tuple of ints

          * With default options, counts contains:                  (nonsyn, syn)
          * If cpg is True (but nonsense is not), counts contains:  (nonsyn, syn, cpg_nonsyn, cpg_syn)
          * If nonsense is True (but cpg is not), counts contains:  (nonsyn, syn, nonsense)
          * If both cpg and nonsense are true, counts contains:     (nonsyn, syn, nonsense, cpg_nonsyn, cpg_syn, cpg_nonsense)
    """
    cds=cds.upper()
    cds2=cds2.upper()  

    syn=0
    nonsyn=0
    nonsensen=0
    if cpg:
        cpg_syn=0
        cpg_nonsyn=0
        cpg_nonsense=0

    dna_letts=set('ACTG')

    if len(cds)%3!=0 or len(cds2)%3!=0:
        raise Exception("count_coding_changes ERROR the sequences must be composed of codons (length multiple of 3)")
    if len(cds)!=len(cds2):
        raise Exception("count_coding_changes ERROR the sequences do not have the same length")

    for i_codon in range(len(cds)//3):
        ## cycling codons
        codon=cds[i_codon*3:i_codon*3+3]
        if not all([lett in dna_letts for lett in codon] ):       
            if not silent:
                warnings.warn(f'count_coding_changes WARNING skipping cds codon n.{i_codon+1} : {codon}')
            continue
        codon2=cds2[i_codon*3:i_codon*3+3]
        if not all([lett in dna_letts for lett in codon2] ):       
            if not silent:
                warnings.warn(f'count_coding_changes WARNING skipping cds2 codon n.{i_codon+1} : {codon}')
            continue

        if ( ('-' in codon and codon!='---') or
             ('-' in codon2 and codon2!='---') ):
             raise Exception("count_coding_changes ERROR the sequences must be aligned by codon, i.e. the gaps must be in groups of three")
         
        if '-' in codon or '-' in codon2: 
            #skipping gap position
            continue

        if codon!=codon2:

            if translate(codon)!=translate(codon2):          
                if "*" in (translate(codon), translate(codon2)):
                    is_nonsyn=2 #nonsense mutation
                else:
                    is_nonsyn=1
            else:
                is_nonsyn=0
                
            for i_within_codon in range(3):
                if codon[i_within_codon]!=codon2[i_within_codon]:
                    i_cds=i_codon*3+i_within_codon
                    nt=cds[i_cds]

                    if not is_nonsyn:    syn+=1
                    elif is_nonsyn==2:   nonsensen+=1
                    else:                nonsyn+=1

                    if cpg:
                        is_cpg= (   nt == 'G' and (  (  i_cds+1<len(cds)  and   cds[i_cds+1] =='C' )  or (i_cds!=0 and cds[i_cds-1] =='C' )    )     ) or \
                            (   nt == 'C' and (  (  i_cds+1<len(cds)  and   cds[i_cds+1] =='G' )  or (i_cds!=0 and cds[i_cds-1] =='G' )    )     )        # # G and the next or previous is C OR #C and the
                    
                        if is_cpg: 
                            if not is_nonsyn:    cpg_syn+=1
                            elif is_nonsyn==2:   cpg_nonsense+=1
                            else:                cpg_nonsyn+=1                
    
    if cpg:
        if nonsense:
            return ( nonsyn, syn, nonsensen, cpg_nonsyn,  cpg_syn, cpg_nonsense )
        else:
            return ( (nonsyn+nonsensen), syn, (cpg_nonsyn+cpg_nonsense),  cpg_syn )
    elif nonsense:
        return ( nonsyn, syn, nonsensen )
    else:
        ## default return value
        return ( (nonsyn+nonsensen), syn )

def count_unique_changes(cds, other_cds_list, silent=True, cpg=False, nonsense=False):
    """ Counts the number of unique synonymous and non-synonymous changes between one reference CDS sequence and a set of others.

     This function is instrumental for any dN/dS analysis.
     Sequences should be aligned using gaps "-". Gapped positions are not counted.
     Codons with any character different from ACTG (for example, N) are skipped and a message is printed to stderr.
     If the same mutation is observed between cds and multiple sequences in other_cds_list, it is counted only once.

    silent : bool, defaults to True
        do not print warning messages for non-ATGC characters

    cpg : bool
        if True, separate counts for CpG sites only are provided, altering the structure of the returned tuple (see below).
        Note that only cds is inspected to identity CpG sites, so with cpg=True the output may be different for (cds, cds2) and (cds2, cds)
        Note that the total counts are still provided. To obtain the number of nonCpG changes, subtract the number of CpG sites from the total number

    nonsense : bool
        if True, non-sense mutations (i.e. stop codons or mutations to a stop codon) are not counted as non-synonymous, and are returned separately, 
        altering the structure of the returned tuple (see below).

    Returns
    -------
    counts : tuple of ints

        * With default options, counts contains:                  (nonsyn, syn)
        * If cpg is True (but nonsense is not), counts contains:  (nonsyn, syn, cpg_nonsyn, cpg_syn)
        * If nonsense is True (but cpg is not), counts contains:  (nonsyn, syn, nonsense)
        * If both cpg and nonsense are true, counts contains:     (nonsyn, syn, nonsense, cpg_nonsyn, cpg_syn, cpg_nonsense)
    """
    cds=cds.upper()
    position_to_change={}  ### hash keeping track of changes we already saw   k: position (0based) -> list of other nts observed in any of other_cds

    syn=0
    nonsyn=0
    nonsensen=0
    if cpg:
        cpg_syn=0
        cpg_nonsyn=0
        cpg_nonsense=0

    dna_letts=set('ACTG')
    for cds2 in other_cds_list:
        cds2=cds2.upper()  

        if len(cds)%3!=0 or len(cds2)%3!=0:
            raise Exception("count_coding_changes ERROR the sequences must be composed of codons (length multiple of 3)")
        if len(cds)!=len(cds2):
            raise Exception("count_coding_changes ERROR the sequences do not have the same length")

        for i_codon in range(len(cds)//3):
            ## cycling codons
            codon=cds[i_codon*3:i_codon*3+3]
            if not all([lett in dna_letts for lett in codon] ):       
                if not silent:
                    warnings.warn(f'count_coding_changes WARNING skipping cds codon n.{i_codon+1} : {codon}')
                continue
            codon2=cds2[i_codon*3:i_codon*3+3]
            if not all([lett in dna_letts for lett in codon2] ):       
                if not silent:
                    warnings.warn(f'count_coding_changes WARNING skipping cds2 codon n.{i_codon+1} : {codon}')
                continue

            if ( ('-' in codon and codon!='---') or
                 ('-' in codon2 and codon2!='---') ):
                 raise Exception("count_coding_changes ERROR the sequences must be aligned by codon, i.e. the gaps must be in groups of three")

            if '-' in codon or '-' in codon2: 
                #skipping gap position
                continue

            if codon!=codon2:

                if translate(codon)!=translate(codon2):          
                    if "*" in (translate(codon), translate(codon2)):
                        is_nonsyn=2 #nonsense mutation
                    else:
                        is_nonsyn=1
                else:
                    is_nonsyn=0

                for i_within_codon in range(3):
                    if codon[i_within_codon]!=codon2[i_within_codon]:
                        i_cds=i_codon*3+i_within_codon
                        nt=cds[i_cds]
                        nt2=cds2[i_cds]
                        # determining if the change is uniq
                        if i_cds in position_to_change and nt2 in position_to_change[i_cds]:
                            continue  #skipping non-uniq
                        # it is uniq
                        if i_cds not in position_to_change:
                            position_to_change[i_cds]=[]
                        position_to_change[i_cds].append(nt2)
                        
                        if not is_nonsyn:    syn+=1
                        elif is_nonsyn==2:   nonsensen+=1
                        else:                nonsyn+=1

                        if cpg:
                            is_cpg= (   nt == 'G' and (  (  i_cds+1<len(cds)  and   cds[i_cds+1] =='C' )  or (i_cds!=0 and cds[i_cds-1] =='C' )    )     ) or \
                                (   nt == 'C' and (  (  i_cds+1<len(cds)  and   cds[i_cds+1] =='G' )  or (i_cds!=0 and cds[i_cds-1] =='G' )    )     )        # # G and the next or previous is C OR #C and the

                            if is_cpg: 
                                if not is_nonsyn:    cpg_syn+=1
                                elif is_nonsyn==2:   cpg_nonsense+=1
                                else:                cpg_nonsyn+=1                

    if cpg:
        if nonsense:
            return ( nonsyn, syn, nonsensen, cpg_nonsyn,  cpg_syn, cpg_nonsense )
        else:
            return ( (nonsyn+nonsensen), syn, (cpg_nonsyn+cpg_nonsense),  cpg_syn )
    elif nonsense:
        return ( nonsyn, syn, nonsensen )
    else:
        ## default return value
        return ( (nonsyn+nonsensen), syn )
    




