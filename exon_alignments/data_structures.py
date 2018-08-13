"""
    The data structures used to model the transcripts

    Supported for both fasta flies and local sqlite db
"""


import sqlite3
from Bio import SeqIO
from tqdm import tqdm
import os


class Transcript_From_Fasta(object):
    """
        Data Structure for Transcript

        Works with fasta files
    """

    def __init__(self, fasta_file, speices_name=None, gene_id=None, gene_name=None, transcript_id=None, transcript_name=None):
        """
        Constructor function

            Loads the fasta_file
            and extracts the exon sequences

            Also contains the metadata (optional)
        
        Arguments:
            fasta_file {str} -- path to fasta file
        
        Keyword Arguments:
            speices_name {str} -- species name (default: {None})
            gene_id {str} -- ensembl id of the corresponding gene (default: {None})
            gene_name {str} -- name of the corresponding gene (default: {None})
            transcript_id {str} -- ensembl id of the transcript (default: {None})
            transcript_name {str} -- name of the transcript (default: {None})
        """
                                            # loading all the data from the arguments
        self.speices_name = speices_name
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.transcript_id = transcript_id
        self.transcript_name = fasta_file.split('_')[-3] + fasta_file.split('_')[-2]  # Extract the transcript name 
                                                                                      # from name of the fasta file

        self.exon_sequence = []         # lists to store exon sequences and ids
        self.exon_id = []
        
        exons = SeqIO.parse(fasta_file, "fasta")  # parse all the exons from the fasta file

        for exon in exons:                      # storing the data in the lists
            self.exon_id.append(exon.description.split(' ')[1])
            self.exon_sequence.append(exon.seq)
            

class Transcript_From_DB(object):
    """
        Data structure for Transcript 

        Works with sqlite database
    """


    def __init__(self, path_to_db, transcript_id, speices_name=None, gene_id=None, gene_name=None,include_UTR=False):
        """
        Constructor Function

            Connects to the database
            and extracts the exon sequences

            Also contains the metadata (optional)
        
        Arguments:
            path_to_db {str} -- path to the database
            transcript_id {str} -- ensembl id of the transcript
        
        Keyword Arguments:
            speices_name {str} -- name of the species (default: {None})
            gene_id {str} -- ensembl id of the corresponding gene (default: {None})
            gene_name {str} -- name of the corresponding gene (default: {None})
            include_UTR {bool} -- include UTR with the exon ? (default: {False})
        """

                                            # loading all the data from the arguments
        self.speices_name = speices_name
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.transcript_id = transcript_id

        
        conn = sqlite3.connect(path_to_db)  # connecting to the database
        c = conn.cursor()

        # Extract all the necessary information related to the transcript from the database
        c.execute('SELECT Exon_ID, seq, Exon_No, start_phase, end_phase FROM Transcripts'+transcript_id[-3:]+" WHERE transcript_id='"+ transcript_id +"'")
        all_rows = c.fetchall()
                                            # initialize lists to store the required information
        start_phases = [0] * len(all_rows)
        end_phases = [0] * len(all_rows)
        self.exon_sequence = [0] * len(all_rows)
        self.exon_id = [0] * len(all_rows)
        self.utr_sequence = [0] * len(all_rows)
        self.utr_pos = [0] * len(all_rows)

                                # iterating through all the exons 
                                # here a point to note is that when we load the exons from the db
                                # it may be the case that they are not sorted even if we do ORDER BY
                                # as all the data in db is in string and 10 is lexicographically smaller than 2 
                                # thus we considered the exon number to keep the exons sorted
                                
        for row in all_rows:

            self.exon_id[int(row[2])-1] = row[0]
            self.exon_sequence[int(row[2])-1] = row[1]
            start_phases[int(row[2])-1] = int(row[3])
            end_phases[int(row[2])-1] = int(row[4])

        if(include_UTR):                # if UTR is not needed not to be filtered, we return the Transcript object as it is

            return
        
        total_transcript_len = 0            # used to identify the coding region from the exons
        
                                            # useful information for UTR removal
        c.execute("SELECT coding_region_start, coding_region_end FROM GeneTranscripts WHERE transcript_id='"+ transcript_id +"'")
        all_rows = c.fetchall()

        coding_region_start = int(all_rows[0][0])
        coding_region_end = int(all_rows[0][1])
                                                        # iterating thorugh all the exons
        for exon_no in range(len(self.exon_sequence)):

            total_transcript_len += len(self.exon_sequence[exon_no])         # adding the sequence len to the counter

            if(start_phases[exon_no] == -1 and end_phases[exon_no] == -1):   # non coding, the entire seq is UTR
                
                self.utr_sequence[exon_no] = self.exon_sequence[exon_no]
                self.utr_pos[exon_no] = -1
                self.exon_sequence[exon_no] = ''

            elif(start_phases[exon_no] == -1 and end_phases[exon_no] != -1):       # the beginning has some UTR 

                left_UTR_len = total_transcript_len - coding_region_start           # finding the UTR region 
                                                                                                                        # sequence splicing
                self.utr_sequence[exon_no] = self.exon_sequence[exon_no][:len(self.exon_sequence[exon_no]) - left_UTR_len - 1] 
                self.utr_pos[exon_no] = -1          # for the visualization tool, -1 means UTR on the left
                self.exon_sequence[exon_no] = self.exon_sequence[exon_no][len(self.exon_sequence[exon_no]) - left_UTR_len - 1 :]

            elif(start_phases[exon_no] != -1 and end_phases[exon_no] == -1):     # the ending has some UTR 

                right_UTR_len = total_transcript_len - coding_region_end            # finding the UTR region 
                                                                                                                        # sequence splicing
                self.utr_sequence[exon_no] = self.exon_sequence[exon_no][len(self.exon_sequence[exon_no]) - right_UTR_len :]
                self.utr_pos[exon_no] = 1          # for the visualization tool, 1 means UTR on the right
                self.exon_sequence[exon_no] = self.exon_sequence[exon_no][  : len(self.exon_sequence[exon_no]) - right_UTR_len ]

        conn.close()            # closing the connection to the database


def main():
    
    # some example codes 

    conn = sqlite3.connect(os.path.join('..','data_acquistion','data','db.db'))
    c = conn.cursor()

    c.execute("SELECT Transcript_ID FROM GeneTranscripts WHERE biotype='protein_coding' AND gencode_basic='1'")
    all_rows = c.fetchall()

    maxx = 0
    maxx_exon_id = ''
    maxx_transcript_id = ''

    for row in tqdm(all_rows):

        transcript = Transcript_From_DB(os.path.join('..','data_acquistion','data','db.db'),row[0])

        for i in range(len(transcript.exon_sequence)):

            if maxx < len(transcript.exon_sequence[i]):

                maxx = len(transcript.exon_sequence[i])
                maxx_exon_id = transcript.exon_id[i]
                maxx_transcript_id = row[0]

        tqdm.write(str(maxx) + ' '+ maxx_exon_id + ' ' + maxx_transcript_id)

    print(maxx)
    

if __name__ == '__main__':
    main()

    