"""
    This is a very basic initial implementation 
    of exon similarity idea. It is currently more for
    understanding the patterns of the data than 
    for the actual algorithm.

    This is a very simplified version, which only works
    for equal number of exons

    Warning : Since now we have the dynamic programming 
    implementation, this simple version may not be compatible 
    with the flask server. But may be used independently 
"""


from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import os

class Transcript(object):

    """
        Data structure for a transcript

        For now it's just a simple list 
        of exon sequences
    """

    def __init__(self, fasta_file, speices_name='None', gene_id='None', gene_name='None', transcript_id='None', transcript_name='None'):

        """
            Loads the fasta_file
            and extracts the exon sequences

            Also contains the metadata (optional)

            :param fasta_file: string , path to the fasta file
            :param speices_name: string, name of the species, optional (default value : None)
            :param gene_id: string, id of the gene, optional (default value : None)
            :param gene_name: string, name of the gene, optional (default value : None)
            :param transcript_id:string, id of the transcript, optional (default value : None)
            :param transcript_name: string, name of the transcript, optional (default value : None)
        """

        self.speices_name = speices_name
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.transcript_id = transcript_id
        self.transcript_name = fasta_file.split('_')[-3] + fasta_file.split('_')[-2]  # getting the transcript name from the file name

        self.exon_sequence = []
        self.exon_id = []
        
        # parsing the fasta file
        exons = SeqIO.parse(fasta_file, "fasta")
        for exon in exons:
            self.exon_id.append(exon.description.split(' ')[1])
            self.exon_sequence.append(exon.seq)
            

def exon_similarity(exon1_seq, exon2_seq, match_score = 1, mismatch_penalty = 0, gap_open = 0, gap_extend = 0):
    
    """
        Computes the exon similarity using global alignment

        :param exon1_seq: string, the first exon sequnce 
        :param exon2_seq: string, the second exon sequnce 
        :param match_score: float, score for matching, optional (default 1)
        :param mismatch_penalty: float, penalty for mismatch, optional (default 0)
        :param gap_open:  float, penalty for opening a gap, optional (default 0)
        :param gap_extend:  float, penalty for continuing a gap, optional (default 0)

        :return alignment_score: float, normalized alignment score, up to two decimal points
    """

    alignment_score = round(pairwise2.align.globalms(exon1_seq, exon2_seq, match_score, mismatch_penalty, gap_open, gap_extend, score_only=True, one_alignment_only=True) / max(len(exon1_seq), len(exon2_seq)) , 2 )          

    return alignment_score

def process_transcript(query_transcript, transcripts):

    """
        Computes similarity of two transcripts, by assuming they have equal number of exons
        
        :param query_transcript: Transcript object, query transcript
        :param transcripts: List of Transcript objects, ortholog transcripts

        :return out_str: string, the output string
    """
    
    out_str = ''  # the output string 
    out_str += query_transcript.transcript_id+'\t'+str(len(query_transcript.exon_sequence))+' Exons'+'\n'

    # heading
    out_str += 'Exon Similarity Scores\tTranscript Score\tMouse Transcript\tNumber of exons'+'\n'
    scores = []
    
    for transcript in transcripts:

        li = []           # list for storing the similarity scores
        total = 0         # combined similarity score  

        for i in range( min ( len(query_transcript.exon_sequence), len(transcript.exon_sequence) ) ):
            alignment_score = exon_similarity(query_transcript.exon_sequence[i], transcript.exon_sequence[i])  # get the similarity measure
            li.append(alignment_score)  # add to the list
            total += alignment_score    # combine the scores

        scores.append(li)
        
        total /= len(query_transcript.exon_sequence)  # normalize by the total number of exons in query sequence

        # add these information to the output string
        out_str += str(li)+'\t'+str(round(total,2))+'\t'+ transcript.transcript_id+'\t'+str(len(transcript.exon_sequence))+'\n'

    return out_str



def main():
    
    '''
        Sample execution of the script

        It takes fasta files as input
    '''

    # specify the path for human transcript
    human_transcript_path = "./exon_fasta_files/NRXN2/human/Homo_sapiens_NRXN2_209_sequence.fa"
    human_transcript = Transcript(human_transcript_path)
    
    # load all the mouse transcripts
    mouse_transcript_paths = next(os.walk('exon_fasta_files/NRXN2/mouse'))[2]
    
    mouse_transcripts = [Transcript("exon_fasta_files/NRXN2/mouse/"+fasta_file) for fasta_file in mouse_transcript_paths]
     

    # process the transcripts
    print(process_transcript(human_transcript,mouse_transcripts))


if __name__ == '__main__':
    main()