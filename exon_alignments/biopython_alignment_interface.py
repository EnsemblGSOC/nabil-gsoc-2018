"""
    We'll use the BioPython alignment tools if MUSCLE
    fails or causes an error
"""


from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqRecord import SeqRecord
from io import StringIO
from Bio import AlignIO
from Bio import pairwise2
import sys
import os
from blosum_matrices import get_raw_blosum_matrix


def protein_similarity_wrapper_biopython(exon1_seq, exon1_id, exon2_seq, exon2_id,  blosum_mat="blosum62", translation_table="Standard"):
    """
    BioPython Alignment interface for Protein similarity
    
    Arguments:
        exon1_seq {str} -- sequence of 1st exon
        exon1_id {str} -- id of 1st exon
        exon2_seq {str} -- sequence of 2nd exon
        exon2_id {str} -- id of 2nd exon
    
    Keyword Arguments:
        blosum_mat {str} -- name of the blosum matrix (default: {"blosum62"})
        translation_table {str} -- DNA translation table to follow (default: {"Standard"})
    
    Returns:
        list -- alignments obtained from BioPython  
    """
    
    # convert to biopython sequence
    coding_dna_1 = Seq(exon1_seq, generic_dna)
    coding_dna_2 = Seq(exon2_seq, generic_dna)
                                                    
    # translate the dna to protein
    translated_protein_1 = coding_dna_1.translate(table=translation_table, to_stop=True)
    translated_protein_2 = coding_dna_2.translate(table=translation_table, to_stop=True)
                                                    
    # converting the sequence to a string
    protein1_record = str(translated_protein_1)
    protein2_record = str(translated_protein_2)

    blosum_matrix = get_raw_blosum_matrix(blosum_mat)

    alignments = pairwise2.align.globaldx(protein1_record, protein2_record, blosum_matrix)

    return [alignments[0][0], alignments[0][1]]


def transcript_similarity_wrapper_biopython(exon1_seq, exon1_id, exon2_seq, exon2_id, match_score, mismatch_penalty, gap_start, gap_extend):
    """
    BioPython Alignment interface for Transcript similarity

    Arguments:
        exon1_seq {str} -- sequence of 1st exon
        exon1_id {str} -- id of 1st exon
        exon2_seq {str} -- sequence of 2nd exon
        exon2_id {str} -- id of 2nd exon
    
    Returns:
        list -- alignments obtained from BioPython 
    """

    alignments = pairwise2.align.globalms(exon1_seq, exon2_seq, match_score, mismatch_penalty, gap_start, gap_extend)

    return [alignments[0][0], alignments[0][1]]


def main():
    
    pass


if __name__ == '__main__':
    main()
