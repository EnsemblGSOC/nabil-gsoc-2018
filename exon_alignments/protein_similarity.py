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

blosums = { "blosum30" : blosum30,
            "blosum35" : blosum35,
            "blosum40" : blosum40,
            "blosum45" : blosum45,
            "blosum50" : blosum50,
            "blosum55" : blosum55,
            "blosum60" : blosum60,
            "blosum62" : blosum62,
            "blosum65" : blosum65,
            "blosum70" : blosum70,
            "blosum75" : blosum75,
            "blosum80" : blosum80,
            "blosum85" : blosum85,
            "blosum90" : blosum90,
            "blosum95" : blosum95,
            "blosum100" : blosum100 }

muscle_exe_path = "~/muscle3.8.31_i86linux64"

def muscle():

    in_file = "../exon_alignments/input.fasta"
    out_file = "../exon_alignments/output.fasta"
    muscle_cline = MuscleCommandline(muscle_exe_path, input=in_file, out=out_file, clw=True)
    muscle_cline()

    align = AlignIO.read(out_file, "clustal")

    out = []

    for record in align:
        out.append(str(record.seq))
        print(str(record.seq))
    

    return out

def protein_similarity_wrapper(exon1_seq, exon1_id, exon2_seq, exon2_id, translation_table="Vertebrate Mitochondrial",blosumMatrix="blosum62", gap_open=-1, gap_extend=-0.5):

    coding_dna_1 = Seq(exon1_seq, generic_dna)
    coding_dna_2 = Seq(exon2_seq, generic_dna)

    translated_protein_1 = coding_dna_1.translate(table=translation_table, to_stop=True)
    translated_protein_2 = coding_dna_2.translate(table=translation_table, to_stop=True)

    protein1_record = SeqRecord(translated_protein_1, id=exon1_id, description="")
    protein2_record = SeqRecord(translated_protein_2, id=exon2_id, description="")

    seq_records = [protein1_record, protein2_record]

    SeqIO.write(seq_records, "../exon_alignments/input.fasta", "fasta")

    out = muscle()

    return out


    #alignments = pairwise2.align.globalds(tranlated_protein_1, tranlated_protein_2, blosums[blosumMatrix], gap_open, gap_extend, one_alignment_only=True)

    print(len(alignments))
    #out = (pairwise2.format_alignment(*alignments[0]))
    
    #return out
    





def main():
        
    
    protein_similarity(exon1,exon2)

if __name__ == '__main__':
    main()