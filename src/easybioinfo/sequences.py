from Bio.Data import CodonTable
from Bio.Seq import Seq, MutableSeq
from Bio import SeqIO, AlignIO

__all__=['codon_table_standard_with_sec', 'Seq', 'MutableSeq', 'CodonTable', 'SeqIO']

#### codon_table_standard_with_sec:
## to be used with BioSeqObject.translate(table= codon_table_standard_with_sec)
# to translate UGA to Sec
codon_table_standard_with_sec = CodonTable.unambiguous_dna_by_name["Standard"]
codon_table_standard_with_sec.forward_table['TGA']='U'
codon_table_standard_with_sec.stop_codons=['TAA', 'TAG']
codon_table_standard_with_sec.protein_alphabet+='U'
