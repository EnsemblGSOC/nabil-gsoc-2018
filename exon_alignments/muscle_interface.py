"""
    MUSCLE interface
"""


from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from io import StringIO
from Bio import AlignIO
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import *
import sys
import os

fp = open('muscle_path.txt','r')
muscle_exe_path = fp.read()
fp.close()


def muscle():
    """
    Interface for MUSCLE alignment
    
    reads the sequence from 'input.fasta'
    writing the resulting alignment in 'output.fasta'

    Returns:
        list -- alignments obtained from MUSCLE 
    """


    in_file = os.path.join("..","exon_alignments","input.fasta")
    out_file = os.path.join("..","exon_alignments","output.fasta")

    # muscle command line interface, CCLUSTAL W mode
    muscle_cline = MuscleCommandline(muscle_exe_path, input=in_file, out=out_file, clw=True)
    muscle_cline()

    align = AlignIO.read(out_file, "clustal")

    out = []

    for record in align:
        out.append(str(record.seq))
        
    return out


def protein_similarity_wrapper(exon1_seq, exon1_id, exon2_seq, exon2_id, translation_table="Standard"):
    """
    MUSCLE interface for Protein similarity
    
    Arguments:
        exon1_seq {str} -- sequence of 1st exon
        exon1_id {str} -- id of 1st exon
        exon2_seq {str} -- sequence of 2nd exon
        exon2_id {str} -- id of 2nd exon
    
    Keyword Arguments:
        translation_table {str} -- DNA translation table to follow (default: {"Standard"})
    
    Returns:
        list -- alignments obtained from MUSCLE 
    """
    
    # convert to biopython sequence
    coding_dna_1 = Seq(exon1_seq, generic_dna)
    coding_dna_2 = Seq(exon2_seq, generic_dna)
                                                    
    # translate the dna to protein
    translated_protein_1 = coding_dna_1.translate(table=translation_table, to_stop=True)
    translated_protein_2 = coding_dna_2.translate(table=translation_table, to_stop=True)
                                                    
    # converting the sequence to a record object for writing a fasta file
    protein1_record = SeqRecord(translated_protein_1, id=exon1_id, description="")
    protein2_record = SeqRecord(translated_protein_2, id=exon2_id, description="")

    seq_records = [protein1_record, protein2_record]    

    # writing the fasta file, to be used by MUSCLE
    SeqIO.write(seq_records, os.path.join("..","exon_alignments","input.fasta"), "fasta")

    out = muscle()

    return out


def transcript_similarity_wrapper(exon1_seq, exon1_id, exon2_seq, exon2_id):
    """
    MUSCLE interface for Transcript similarity

    Arguments:
        exon1_seq {str} -- sequence of 1st exon
        exon1_id {str} -- id of 1st exon
        exon2_seq {str} -- sequence of 2nd exon
        exon2_id {str} -- id of 2nd exon
    
    Returns:
        list -- alignments obtained from MUSCLE 
    """

    # convert to biopython sequence
    coding_dna_1 = Seq(exon1_seq, generic_dna)
    coding_dna_2 = Seq(exon2_seq, generic_dna)

    # converting the sequence to a record object for writing a fasta file
    transcript1_record = SeqRecord(coding_dna_1, id=exon1_id, description="")
    transcript2_record = SeqRecord(coding_dna_2, id=exon2_id, description="")

    seq_records = [transcript1_record, transcript2_record]

    # writing the fasta file, to be used by MUSCLE
    SeqIO.write(seq_records, os.path.join("..","exon_alignments","input.fasta"), "fasta")

    out = muscle()

    return out


def main():
    
    pass


if __name__ == '__main__':
    main()