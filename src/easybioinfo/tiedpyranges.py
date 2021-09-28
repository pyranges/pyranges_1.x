import os, pyranges, pandas,  pyfaidx

__all__=['get_faidx_handler', 'pyranges']

########################################################
############ General methods

def get_faidx_handler(fastafile):
    """ Returns a pyfaidx.Fasta object to fetch sequences from a fasta file

    Parameters
    ----------
    fastafile : str
        filename or path to fasta fileStart of subregion, 0-based and included. Use a negative int to count from the 3'  (e.g. -1 is the last nucleotide)

    Returns
    -------
    pf : pyfaidx.Fasta
        sequence I/O handler to fetch sequence

    Note
    ----
    This function will create a .fai file if not present. In doing so, it creates a 
     temporary INFILE.tmplnk.fai file first, and upon completion it is moved to INFILE.fai
    
    """
    indexfile=fastafile+'.fai'
    if os.path.isfile(indexfile):
        return pyfaidx.Fasta(fastafile)
    else:
        ## Creating temporary link to fasta, so that
        # index will be created in temporary location
        temp_link= fastafile+'.tmplnk'
        temp_index=fastafile+'.tmplnk.fai'
        for i in [temp_link, temp_index]:
            if os.path.isfile(i):
                os.remove(i)
        os.symlink(  os.path.basename(fastafile),  temp_link   )

        # building index
        service(f'Indexing {fastafile} with faidx...')
        pf=pyfaidx.Fasta(temp_link)
        flush_service()
        
        # done, cleaning up 
        os.rename(temp_index, indexfile)
        os.remove(temp_link)
        pf.filename=fastafile
        pf.faidx.filename=fastafile
        pf.faidx.indexname=indexfile
        return pf

########################################################
############ Methods added to PyRanges:

def cut(self, start=0, length=None, end=None, by=None, is_sorted=False):
    """ Cuts a portion of the intervals, corresponding to genomic subregions of the input (self)

    Parameters
    ----------
    start : int
        Start of subregion, 0-based and included. Use a negative int to count from the 3'  (e.g. -1 is the last nucleotide)

    length : int, default None
        Length of subregion. If None, everything after start is returned

    end : int, default None
        End of subregion. Alternative method to define subregion to providing length

    by : list of str, default None
        intervals are grouped by this/these ID column(s) beforehand, e.g. exons belonging to same transcripts

    is_sorted : bool, default False
        if self is already sorted by Start, set this to True to speed up computation

    Returns
    -------
    PyRanges
        Subregion of self, cut as specified by arguments

    Note
    ----
    If the request goes out of bounds (e.g. requesting 100 nts for a 90nt region), only the existing portion is returned


    Examples
    --------
    >>>p  = pr.from_dict({"Chromosome": [1, 1, 2, 2, 3], 
    ...                   "Strand": ["+", "+", "-", "-", "+"],
    ...                   "Start": [1, 40, 2, 30,  140], 
    ...                   "End": [20, 60, 13, 45, 155], 
    ...                   "transcript_id":["t1", "t1", "t2", "t2", "t3"] })
    >>>p
    +--------------+--------------+-----------+-----------+-----------------+
    |   Chromosome | Strand       |     Start |       End | transcript_id   |
    |   (category) | (category)   |   (int32) |   (int32) | (object)        |
    |--------------+--------------+-----------+-----------+-----------------|
    |            1 | +            |         1 |        20 | t1              |
    |            1 | +            |        40 |        60 | t1              |
    |            2 | -            |         2 |        13 | t2              |
    |            2 | -            |        30 |        45 | t2              |
    |            3 | +            |       140 |       155 | t3              |
    +--------------+--------------+-----------+-----------+-----------------+
    Stranded PyRanges object has 5 rows and 5 columns from 3 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    Get the first 10 nucleotides (at the 5') of *each interval* (each line of the dataframe):
    >>>p.cut(0, 10)
    +--------------+--------------+-----------+-----------+-----------------+
    |   Chromosome | Strand       |     Start |       End | transcript_id   |
    |   (category) | (category)   |   (int32) |   (int32) | (object)        |
    |--------------+--------------+-----------+-----------+-----------------|
    |            1 | +            |         1 |        20 | t1              |
    |            1 | +            |        40 |        60 | t1              |
    |            2 | -            |         2 |        13 | t2              |
    |            2 | -            |        30 |        45 | t2              |
    |            3 | +            |       140 |       155 | t3              |
    +--------------+--------------+-----------+-----------+-----------------+
    Stranded PyRanges object has 5 rows and 5 columns from 3 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    Get the first 10 nucleotides of *each transcript*, grouping exons by transcript_id:
    >>>p.cut(0, 10, by='transcript_id')
    +--------------+--------------+-----------+-----------+-----------------+
    |   Chromosome | Strand       |     Start |       End | transcript_id   |
    |   (category) | (category)   |   (int32) |   (int32) | (object)        |
    |--------------+--------------+-----------+-----------+-----------------|
    |            1 | +            |         1 |        11 | t1              |
    |            2 | -            |        35 |        45 | t2              |
    |            3 | +            |       140 |       150 | t3              |
    +--------------+--------------+-----------+-----------+-----------------+
    Stranded PyRanges object has 3 rows and 5 columns from 3 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    Get the last 20 nucleotides of each transcript:
    >>>p.cut(-20,  by='transcript_id')
    +--------------+--------------+-----------+-----------+-----------------+
    |   Chromosome | Strand       |     Start |       End | transcript_id   |
    |   (category) | (category)   |   (int32) |   (int32) | (object)        |
    |--------------+--------------+-----------+-----------+-----------------|
    |            1 | +            |        40 |        60 | t1              |
    |            2 | -            |        30 |        39 | t2              |
    |            2 | -            |         2 |        13 | t2              |
    |            3 | +            |       140 |       150 | t3              |
    +--------------+--------------+-----------+-----------+-----------------+
    Stranded PyRanges object has 4 rows and 5 columns from 3 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Stra

    Get region from 30 to 330 of each transcript, or their existing subportion:
    >>>p.cut(30, 300, by='transcript_id')
    +--------------+--------------+-----------+-----------+-----------------+
    |   Chromosome | Strand       |     Start |       End | transcript_id   |
    |   (category) | (category)   |   (int32) |   (int32) | (object)        |
    |--------------+--------------+-----------+-----------+-----------------|
    |            1 | +            |        51 |        60 | t1              |
    +--------------+--------------+-----------+-----------+-----------------+
    Stranded PyRanges object has 1 rows and 5 columns from 1 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.
    """
    from_end=[False,False]
    if start<0:
        start=-start
        from_end[0]=True
    if not end is None and end<0:
        end=-end
        from_end[1]=True

    z=self.sort()  if not is_sorted   else self

    if by is None:
        return  pyranges.PyRanges(
            pandas.concat(
                [v.groupby(  pandas.Index( range(len(v)))).apply(
                    _cut_pyranges,
                    start=start, length=length, end=end, from_end=from_end)
                 for k,v in z.items()   ] ))
        
 
    else:
        return  pyranges.PyRanges(
            pandas.concat(
                [v.groupby(  by, sort=False ).apply(
                    _cut_pyranges,
                    start=start, length=length, end=end, from_end=from_end)
                 for k,v in z.items()   ]) )

pyranges.PyRanges.cut=cut
###########

def extend(self, up=0, down=0, by=None, is_sorted=False):
    """ Extend a stranded genomic interval upstream or downstream

    Parameters
    ----------
    up : int, default 0
        How much to extend upstream

    down : int, default 0
        How much to extend downstream

    by : list of str, default None
        intervals are grouped by this/these ID column(s) beforehand, e.g. exons belonging to same transcripts

    is_sorted : bool, default False
        if self is already sorted by Start, set this to True to speed up computation

    Returns
    -------
    PyRanges
        Extended by one or both sides as per specifications

    Note
    ----
    The request may go out of bounds, either below 0 or above the chromosome size. It's best to run pyranges.genomicfeatures.genome_bounds after


    Examples
    --------
    >>>p  = pr.from_dict({"Chromosome": [1, 1, 2, 2, 3], 
    ...                   "Strand": ["+", "+", "-", "-", "+"],
    ...                   "Start": [1, 40, 2, 30,  140], 
    ...                   "End": [20, 60, 13, 45, 155], 
    ...                   "transcript_id":["t1", "t1", "t2", "t2", "t3"] })
    >>>p
    +--------------+--------------+-----------+-----------+-----------------+
    |   Chromosome | Strand       |     Start |       End | transcript_id   |
    |   (category) | (category)   |   (int32) |   (int32) | (object)        |
    |--------------+--------------+-----------+-----------+-----------------|
    |            1 | +            |         1 |        20 | t1              |
    |            1 | +            |        40 |        60 | t1              |
    |            2 | -            |         2 |        13 | t2              |
    |            2 | -            |        30 |        45 | t2              |
    |            3 | +            |       140 |       155 | t3              |
    +--------------+--------------+-----------+-----------+-----------------+
    Stranded PyRanges object has 5 rows and 5 columns from 3 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    Extend by 3nt downstream *each interval* (each line of the dataframe):
    >>>p.extend(down=2)
    +--------------+--------------+-----------+-----------+-----------------+
    |   Chromosome | Strand       |     Start |       End | transcript_id   |
    |   (category) | (category)   |   (int32) |   (int32) | (object)        |
    |--------------+--------------+-----------+-----------+-----------------|
    |            1 | +            |         1 |        22 | t1              |
    |            1 | +            |        40 |        62 | t1              |
    |            2 | -            |         0 |        13 | t2              |
    |            2 | -            |        28 |        45 | t2              |
    |            3 | +            |       140 |       157 | t3              |
    +--------------+--------------+-----------+-----------+-----------------+
    Stranded PyRanges object has 5 rows and 5 columns from 3 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    Extend by 10nt upstream and 3nt downstream *each transcript*, grouping exons by transcript_id:
    >>>p.extend(up=10, down=2, by='transcript_id')
    +--------------+--------------+-----------+-----------+-----------------+
    |   Chromosome | Strand       |     Start |       End | transcript_id   |
    |   (category) | (category)   |   (int32) |   (int32) | (object)        |
    |--------------+--------------+-----------+-----------+-----------------|
    |            1 | +            |        -9 |        20 | t1              |
    |            1 | +            |        40 |        62 | t1              |
    |            2 | -            |         0 |        13 | t2              |
    |            2 | -            |        30 |        55 | t2              |
    |            3 | +            |       130 |       157 | t3              |
    +--------------+--------------+-----------+-----------+-----------------+
    Stranded PyRanges object has 5 rows and 5 columns from 3 chromosomes.
    For printing, the PyRanges was sorted on Chromosome and Strand.

    """
    z=self.sort()  if not is_sorted   else self    
    start_col=z.columns.get_loc('Start')
    end_col=  z.columns.get_loc('End')
    
    if by is None:
        return  pyranges.PyRanges(
            pandas.concat(
                [v.groupby(  pandas.Index( range(len(v)))).apply(
                    _extend_pyranges,
                    up=up, down=down, start_col=start_col, end_col=end_col)
                 for k,v in z.items()   ] ))        
 
    else:
        return  pyranges.PyRanges(
            pandas.concat(
                [v.groupby(  by, sort=False ).apply(
                    _extend_pyranges, 
                    up=up, down=down, start_col=start_col, end_col=end_col)                 
                 for k,v in z.items()   ]) )

pyranges.PyRanges.extend=extend
        
def get_sequence(self, by, faidx_handler=None, fastafile=None):
    """ Fetch the sequence corresponding to grouped intervals (exons of genes).

    Parameters
    ----------

    by : list of str
        intervals are grouped by this/these ID column(s) beforehand, e.g. exons belonging to same transcripts

    faidx_handler : faidx.Fasta, default None
        Sequence I/O handler. Get it with: get_faidx_handler(fastafilename)

    fastafile : str, default None
        Alternatively to faidx_handler, provide path to fastafile and the handler will be derived. 
        Note that the file will be indexed if no .fai file is found

    Returns
    -------
    Pandas.DataFrame
        With the column 'Sequence', indexed by all columns in the 'by' argument
    """
    
    if faidx_handler is None:
        faidx_handler=get_faidx_handler(fastafile)
        
    return(pandas.concat(
        [v.groupby( by, sort=False, group_keys=False).apply(
                    _get_seq, faidx_handler)
         for k,v in self.items()]))

pyranges.PyRanges.get_sequence=get_sequence


    

########################################################
############ END of methods added to PyRanges
########################################################

            
def _get_seq(p, faidx_handler):
    if not len(p): return ''
    positions=[]
    revcomp='Strand' in p.columns and p.Strand=='-'
    for r in p.itertuples():
        positions.append([r.Start+1, r.End]  )   #0-based [,)  to 1-based [,]
    seq=faidx_handler.get_spliced_seq( p.Chromosome.iloc[0],
                                       positions,
                                       rc=revcomp ).seq
    #write(seq, how='green')
    return( pandas.Series(seq, index=['Sequence']) )


def _extend_pyranges(p, up=0, down=0, start_col=2, end_col=3):
    # strand= p.Strand.unique()
    # if len(strand)!=1 or not strand[0] in '+-':
    #     raise Exception(f"extend ERROR pyranges must be stranded and unique per group! {strand}" )
    strand=p.Strand.at[0]
    if not strand in '+-':
        raise Exception(f"extend ERROR pyranges must be stranded! {strand}" )
    
    if up:
        if strand=='+':
            p.iat[0, start_col] -= up
        else:
            p.iat[len(p)-1, end_col] += up
            
    if down:
        if strand=='+':
            p.iat[ len(p)-1, end_col] += down
        else:
            p.iat[0, start_col] -= down
            
    return(p)

def _cut_pyranges(p,  start=0,  length=None, end=None, from_end=[False,False]):
    """ Input p is a dataframe from a pyranges object grouped by ID --> a single tied pyrange object
    Start is 0-based included, end is 0-based excluded
    If length is None -> all there is
    If from_end:   start=1 corresponds to python [-1] etc
    """
    strand=p.Strand.at[0]
    if not strand in '+-':
        raise Exception(f"cut ERROR pyranges must be stranded! {strand}" )

    
    if strand=='-':
        # reversing so cumulative length is computed correctly
        p=p[::-1]
    p['_lengs']=p.End-p.Start
    p['_cumlen']= p._lengs.cumsum()
    totlen=p._cumlen.iloc[-1]

    #########
    ## Handling of input arguments
    if length is None and end is None:
        length=totlen
    if from_end[0]:
        start=totlen-start
    if not end is None and from_end[1]:
        end=totlen-end
    if end is None:
        end=start+length
    #########

    if    strand=='+':
        p['_dst']=( p._lengs + start - p._cumlen ).clip(lower=0)
        p['_den']=( end - p._cumlen ).clip(upper=0)
    else:
        p['_dst']=( -(end - p._cumlen) ).clip(lower=0)
        p['_den']=( -(p._lengs + start - p._cumlen ) ).clip(upper=0)

    selector=(p._dst < p._lengs) & (p._lengs > -p._den)

    if selector.empty:
        #raise Exception('ERROR subseq went out of boundaries!')
        return pandas.DataFrame()

    p=p [ selector ]
    p.Start=p.Start + p._dst
    p.End  =p.End   + p._den

    outdf=p.drop( ['_lengs', '_cumlen', '_dst', '_den'], axis=1)
    return outdf
