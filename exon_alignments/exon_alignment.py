'''
    This is a very basic initial implementation 
    of my current idea. It is currently more for
    understanding the patterns of the data than 
    for the actual algorithm.
'''


from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import os

class Transcript(object):

    '''
        Data structure for a transcript

        For now it's just a simple list 
        of exon sequences
    '''

    def __init__(self, fasta_file, speices_name=None, gene_id=None, gene_name=None, transcript_id=None, transcript_name=None):

        '''
            Loads the fasta_file
            and extracts the exon sequences

            Also contains the metadata (optional)
        '''

        self.speices_name = speices_name
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.transcript_id = transcript_id
        self.transcript_name = fasta_file.split('_')[-3] + fasta_file.split('_')[-2]

        self.exon_sequence = []
        self.exon_id = []
        
        exons = SeqIO.parse(fasta_file, "fasta")
        for exon in exons:
            self.exon_id.append(exon.description.split(' ')[1])
            self.exon_sequence.append(exon.seq)
            


def process_transcript(query_transcript, transcripts):

    '''
        Roughly applies a very simplified version
        of the algorithm I proposed.
    '''

    print(query_transcript.transcript_name)

    print('Exon Similarity Scores - Transcript Score - Transcript Name - Number of exons')
    scores = []
    
    for transcript in transcripts:

        li = []
        total = 0

        for i in range( min ( len(query_transcript.exon_sequence), len(transcript.exon_sequence) ) ):
            alignment_score = round(pairwise2.align.localxx(query_transcript.exon_sequence[i], transcript.exon_sequence[i], score_only=True, one_alignment_only=True) / max(len(query_transcript.exon_sequence[i]), len(transcript.exon_sequence[i])),2)
            li.append(alignment_score)
            total += alignment_score

        scores.append(li)
        
        print(li,' - ', round(total,2),' - ', transcript.transcript_name, '-', len(transcript.exon_sequence))




def main():
    
    '''
        Currently I'm just comparing a human transcript with the ortholog
        mouse transcripts.

        Currently it takes fasta files as input, which will soon be 
        migrated to a database.
    '''

    # specify the path for human transcript
    human_transcript_path = "exon_fasta_files/NRXN2/human/Homo_sapiens_NRXN2_201_sequence.fa"
    human_transcript = Transcript(human_transcript_path)
    
    # loads all the mouse transcripts
    mouse_transcript_paths = next(os.walk('exon_fasta_files/NRXN2/mouse'))[2]
    
    mouse_transcripts = [Transcript("exon_fasta_files/NRXN2/mouse/"+fasta_file) for fasta_file in mouse_transcript_paths]
     
    process_transcript(human_transcript,mouse_transcripts)

if __name__ == '__main__':
    main()